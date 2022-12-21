%% Marco Iglesias, Universtity of Nottingham, 2022
function [Storage_kidney,Loss_kidney,Level]=physical(Un,Model,prior)
%% This function transforms the parameters from Un into physical quantities 
% i.e. Storage and Loss modulus. It also provides the level-set function
% that defines the cyst/tumour

mask=Model.mask;
b=prior.max_len;
a=prior.min_len;
N=Model.Grid.N;


Level_lengths=Un{1,2};

Storage_mean=Un{1,3};
Storage_RHS=Un{1,4};
Storage_lengths=Un{1,5};

Loss_mean=Un{1,6};
Loss_RHS=Un{1,7};
Loss_lengths=Un{1,8};

NF=prior.no_fields;
%% Compute Level-set function
for idim=1:Model.Grid.dim
    pri.len{idim}=(1/pi*atan(Level_lengths(idim))+0.5)*(b-a)+a;
    pri.len{idim}=pri.len{idim}*ones(N,1);
end
pri.sigma=prior.level(1).sigma; pri.nu=prior.level(1).nu;
Level=grf2D(Model, pri,Un{1,1});

%% Compute Storage and Loss on each region
Storage=cell(NF,1);
Loss=cell(NF,1);
for nf=1:NF
    Storage{nf}=Storage_mean(nf);
    for idim=1:Model.Grid.dim
        Sto_lengths=Storage_lengths(idim+(nf-1)*Model.Grid.dim)*ones(N,1);
        pri.len{idim}=(1/pi*atan(Sto_lengths)+0.5)*(b-a)+a;
    end
    pri.sigma=prior.Storage(nf).sigma; 
    pri.nu=prior.Storage(nf).nu;    
    Storage{nf}=Storage{nf}+grf2D(Model, pri,Storage_RHS(1+(nf-1)*N:N+(nf-1)*N,1));

    Loss{nf}=Loss_mean(nf);
    for idim=1:Model.Grid.dim
        Lo_lengths=Loss_lengths(idim+(nf-1)*Model.Grid.dim)*ones(N,1);
        pri.len{idim}=(1/pi*atan(Lo_lengths)+0.5)*(b-a)+a;
    end
    pri.sigma=prior.Loss(nf).sigma; 
    pri.nu=prior.Loss(nf).nu;    
    Loss{nf}=Loss{nf}+grf2D(Model, pri,Loss_RHS(1+(nf-1)*N:N+(nf-1)*N,1));    
end


Nx=Model.Grid.Nx;
Ny=Model.Grid.Ny;
bet=prior.threshold;

%% Put properties together using the mask
Storage_temp=reshape(exp(Storage{2}),Nx,Ny).*( mask==1)+reshape(exp(Storage{3}),Nx,Ny).*( mask==2);
Storage_tempB=Storage_temp+(-Storage_temp+reshape(exp(Storage{4}),Nx,Ny)).*(reshape(Level,Nx,Ny)>=bet);
Storage_kidney=reshape(exp(Storage{1}),Nx,Ny).*( mask==0)+Storage_tempB.*( mask~=0);

Loss_temp=reshape(exp(Loss{2}),Nx,Ny).*( mask==1)+reshape(exp(Loss{3}),Nx,Ny).*( mask==2);
Loss_tempB=Loss_temp+(-Loss_temp+reshape(exp(Loss{4}),Nx,Ny)).*(reshape(Level,Nx,Ny)>=bet);
Loss_kidney=reshape(exp(Loss{1}),Nx,Ny).*( mask==0)+Loss_tempB.*( mask~=0);

Loss_kidney=reshape(Loss_kidney,Model.Grid.Nx*Model.Grid.Ny,1);
Storage_kidney=reshape(Storage_kidney,Model.Grid.Nx*Model.Grid.Ny,1);


end




