%% Marco Iglesias, Universtity of Nottingham, 2022
function EKI(Un,Model,prior,N_En, data,inv_sqrt_Gamma)
%% This is the code that runs ensemble Kalman inversion 

M=length(data);
f_mean= @(A) mean(A,2);
Un_m = cellfun(f_mean,Un,'un',0);

%% Initialisations
MAX_iter=120;
t=zeros(MAX_iter,1);
Misfit=zeros(MAX_iter,1);
t(1)=0;
Cond=1;
iter=0;
Output=cell(N_En);
Storage=zeros(Model.Grid.N,N_En);
Loss=zeros(Model.Grid.N,N_En);
Indi=zeros(Model.Grid.N,N_En);

num_work=6;
delete(gcp('nocreate'))
fprintf('Number of slots available: %d\n', num_work);
parpool('local', num_work);
Z=zeros(M,N_En);

%% EKI iterations

while (iter<MAX_iter)
    iter=iter+1;    
    %Prediction Step: Run forward model for each parameter in Un and get
    %Storage, Loss and level-set
    parfor en = 1:N_En
        [Output{en},Storage(:,en),Loss(:,en),Indi(:,en)]=...'
            fwd_model(cell_en(Un,en),Model,prior);
    end
    
    meanSto=mean(Storage,2);
    meanLoss=mean(Loss,2);
    VarSto=var(Storage,0,2);
    VarLoss=var(Loss,0,2);
    prob=mean(Indi,2);
    StoMeans=Un{1,3};
    LossMeans=Un{1,6};
    if (iter==1)
        save('Results_first','Model', 'prior','VarLoss', 'VarSto','prob',...
            'meanSto','meanLoss','StoMeans','LossMeans','-v7.3')
    end
    % Compute weighted data misfits
    for en=1:N_En
        Z(:,en)=inv_sqrt_Gamma.*(data-Output{en});
    end

    Z_m=mean(Z,2);
    Misfit(iter)=norm(Z_m)^2/M;

    if (Cond==0)
        break
    end
    %%get regularisation parameter as in
    % Marco Iglesias and Yuchen Yang 2021 Inverse Problems 37 025008
    alpha=mean(vecnorm(Z).^2)/M;
    
    %Check if we have arrive at approximate posterior
    if (t(iter)+1/alpha>1)
        alpha=1/(1-t(iter));
        Cond=0;
    end
    % Update parameters via Kalman formula
    Un=update(Un,Un_m,alpha,Z,Z_m,N_En);
    %update means
    Un_m = cellfun(f_mean,Un,'un',0);

    t(iter+1)=t(iter)+1/alpha;
    
    disp(['Misfit: [',num2str(Misfit(1:iter)'),']'])
    disp(['t: [',num2str(t(1:iter+1)'),']'])

   
end

save('Results_final','Model', 'prior','VarLoss', 'VarSto','prob',...
            'meanSto','meanLoss','StoMeans','LossMeans','-v7.3')


end
function U_en=cell_en(Un,en)
%% Get an ensemble member en from cell Un
U_en=cell(1,length(Un));
for i=1:length(Un)
    U_en{1,i}=Un{1,i}(:,en);
end
end


