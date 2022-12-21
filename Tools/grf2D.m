%% Marco Iglesias, Universtity of Nottingham, 2022
function U=grf2D(Model, prior,xi)
%%We use the SPDE approach to compute Random fields, 
% Here xi is the Right hand side and prior is a structure that has the
% lengthscales, smoothness parameter nu and amplitude scale. 
% for further details see appendix in
% Marco Iglesias et al 2022 Phys. Med. Biol. 67 235003

lx=prior.len{1};
ly=prior.len{2};
sigma=prior.sigma;
nu=prior.nu;
dim=Model.Grid.dim;

alpha2=sigma^2*2^dim*pi^(dim/2)*gamma(nu+dim/2)/gamma(nu);
Grid=Model.Grid;
    
Lambda_x=lx;
Lambda_y=ly;

hx=1/Grid.Nx;
hy=1/Grid.Ny;

f=alpha2*Lambda_x.*Lambda_y*hx*hy;
sqrt_f=sqrt(alpha2*Lambda_x.*Lambda_y)*hx*hy;
%% Get stiffness matrix
A=get_Mats(Grid,Lambda_x,Lambda_y);

if (nu==1)    
    [U,~,~,~,~]  = pcg(A,sqrt(f).*xi,[],1000);    
elseif (nu==2)
    [R,p]=chol(A);
    y=R\xi;
    U=A\(sqrt_f.*y);    
elseif (nu==3)
    [y,~,~,~,~]  = pcg(A,sqrt(f).*xi,[],1000);    
    [U,~,~,~,~]  = pcg(A,hx*hy.*y,[],1000);    
end
end

function A=get_Mats(Grid,Lambda_x,Lambda_y)

Nx=Grid.Nx; Ny=Grid.Ny;  N=Nx*Ny;
hx=1/Nx; hy=1/Ny; 

lam_x=reshape(Lambda_x,Nx,Nx);
lam_y=reshape(Lambda_y,Nx,Nx);


Lx=lam_x.^2;
Ly=lam_y.^2;

fac=1.2;
lambda_left=fac*lam_x(1,:);
lambda_right=fac*lam_x(Nx,:);
lambda_top=fac*lam_y(:,Nx);
lambda_bottom=fac*lam_y(:,1);


tx=hy/hx; TX=zeros(Nx+1,Ny,1);
ty=hx/hy; TY=zeros(Nx,Ny+1,1);
Average_x=0.5*(Lx(1:Nx-1,:)+Lx(2:Nx,:));
TX(2:Nx,:)=Average_x.*tx;
Average_y=0.5*(Ly(:,1:Ny-1)+Ly(:,2:Ny));
TY(:,2:Ny)=Average_y.*ty;
TX2=TX;
TY2=TY;
% 
TX(Nx+1,1:Nx)=(tx./(1/2+lambda_left/hx)).*Lx(1,:);
TX(1,1:Nx)=(tx./(1/2+lambda_right/hx)).*Lx(Nx,:);
TY(1:Nx,1)=(ty./(1/2+lambda_bottom/hx)).*Ly(:,1);
TY(1:Nx,Nx+1)=(ty./(1/2+lambda_top/hx)).*Ly(:,Nx);


x1=reshape(TX(1:Nx,:),N,1); x2=reshape(TX(2:Nx+1,:),N,1);
y1=reshape(TY(:,1:Ny),N,1); y2=reshape(TY(:,2:Ny+1),N,1);

x12=reshape(TX2(1:Nx,:),N,1); x22=reshape(TX2(2:Nx+1,:),N,1);
y12=reshape(TY2(:,1:Ny),N,1); y22=reshape(TY2(:,2:Ny+1),N,1);

DiagVecs=[-y22,-x22,x1+x2+y1+y2,-x12,-y12];
DiagIndx=[-Nx,-1,0,1,Nx];
A=spdiags(DiagVecs,DiagIndx,N,N);

DiagVecs=hx*hy*ones(N,1);
DiagIndx=0;
C=spdiags(DiagVecs,DiagIndx,N,N);
A=A+C;
end

