%% Marco Iglesias, Universtity of Nottingham, 2022
function Un=Get_prior_fields(Pr,Model,N_En)
%% Produce samples from the prior of the parameters and same them in cell Un
b=Pr.max_len;
a=Pr.min_len;

N=Model.Grid.N;
dim=Model.Grid.dim;
Level_lengths=zeros(dim,N_En);
Storage_RHS=zeros(Pr.no_fields*N,N_En);
Storage_means=zeros(Pr.no_fields,N_En);
Storage_lengths=zeros(dim*Pr.no_fields,N_En);

Loss_RHS=zeros(Pr.no_fields*N,N_En);
Loss_means=zeros(Pr.no_fields,N_En);
Loss_lengths=zeros(dim*Pr.no_fields,N_En);

Level_RHS=randn(N,N_En); % Gaussian white noise for the RHS of the SPDE (

X = lhsdesign(N_En,dim);
% Sample from uniform on lenthscales then transform to used transformed
% variables in EKI (to ensure we keep lengthscales in [a,b])
for idim=1:dim
    lengths=Pr.level.len(idim).mean.lim(1)+(Pr.level.len(idim).mean.lim(2)-Pr.level.len(idim).mean.lim(1))*X(:,idim);
    Level_lengths(idim,:)= tan(pi*( (lengths-a)/(b-a)-0.5) );
end

%Samples means for log Storage/Loss, lengthscales, and RHS for SPDE
for ifield=1:Pr.no_fields
    X = lhsdesign(N_En,dim);
    exp_Sto=exp(Pr.Storage(ifield).mean.lim(1))+(exp(Pr.Storage(ifield).mean.lim(2))-exp(Pr.Storage(ifield).mean.lim(1)))* lhsdesign(N_En,1);
    Storage_means(ifield,:)=log(exp_Sto);
    Storage_RHS(1+(ifield-1)*N:N+(ifield-1)*N,:)=randn(N,N_En);
    for idim=1:dim
        lengths=Pr.Storage(ifield).len(idim).lim(1)+(Pr.Storage(ifield).len(idim).lim(2)-Pr.Storage(ifield).len(idim).lim(1))*X(:,idim);
        Storage_lengths(idim+(ifield-1)*dim,:)=tan(pi*( (lengths-a)/(b-a)-0.5) );
    end

end

for ifield=1:Pr.no_fields
    X = lhsdesign(N_En,dim);
    exp_Loss=exp(Pr.Loss(ifield).mean.lim(1))+(exp(Pr.Loss(ifield).mean.lim(2))-exp(Pr.Loss(ifield).mean.lim(1)))*lhsdesign(N_En,1);
    Loss_means(ifield,:)=log(exp_Loss);
    Loss_RHS(1+(ifield-1)*N:N+(ifield-1)*N,:)=randn(N,N_En);
    for idim=1:dim
        lengths=Pr.Loss(ifield).len(idim).lim(1)+(Pr.Loss(ifield).len(idim).lim(2)-Pr.Loss(ifield).len(idim).lim(1))*X(:,idim);
        Loss_lengths(idim+(ifield-1)*dim,:)=tan(pi*( (lengths-a)/(b-a)-0.5) );
    end

end




Un{1,1}=Level_RHS; %RHS for SPDE of level-set random field
Un{1,2}=Level_lengths; %lengthscales for SPDE of level-set random field
Un{1,3}=Storage_means; %mean (constant) for the log-Storage random field
Un{1,4}=Storage_RHS; %RHS for SPDE of log-Storage random field
Un{1,5}=Storage_lengths;%lengthscales for SPDE of of log-Storage random field
Un{1,6}=Loss_means;%mean (constant) for the log-Losse random field
Un{1,7}=Loss_RHS;%RHS for SPDE of log-Loss random field
Un{1,8}=Loss_lengths;%lengthscales for SPDE of of log-Loss random field









