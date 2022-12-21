%% Marco Iglesias, Universtity of Nottingham, 2022
function out=Solver(Model,G,Loss)
%% This function runs the viscoelastic 2D model used in 
% Marco Iglesias et al 2022 Phys. Med. Biol. 67 235003
% given Storage (G) and Loss modulus using the PDE toolbox. 
% The solution (complex displacements in x and y) are interpolated on the
% measurement grid of M x M points

rho=Model.rho;
omega=Model.omega;
load=[0;Model.dispAmplitude];

mu=G+Loss*1i;
nu=Model.nu;
lambda=2*G*nu/(1-2*nu);
model=Model.FM;
mu_lambda=2*mu+lambda;
X=Model.Grid.X;
Y=Model.Grid.Y;

applyBoundaryCondition(model,'neumann','Edge',[2,1,4],'q',0,'g',0);
applyBoundaryCondition(model,'dirichlet','Edge',3,'r',load,'h',[0,0;0,1]);
specifyCoefficients(model,'m',0,'d',0,'c',@ccoeffunction,'a',-rho*omega^2,'f',[0;0]);
results = solvepde(model); 
delete(model.BoundaryConditions);
delete(model.EquationCoefficients);

querypoints=[Model.Grid.X_meas(:),Model.Grid.Y_meas(:)]';
uintrp = interpolateSolution(results,querypoints,[1,2]);
re_ux=real(uintrp (:,1));
re_uy=real(uintrp (:,2));
im_ux=imag(uintrp (:,1));
im_uy=imag(uintrp (:,2));
out=[re_ux;re_uy;im_ux;im_uy];   


function cmatrix = ccoeffunction(location,state)
    n1 = 3;
    nr = numel(location.x);
    x=location.x;
    y=location.y;
    mu_lambda_interp=interp2(X,Y,mu_lambda,x,y);
    mu_interp=interp2(X,Y,mu,x,y);
    lambda_interp=interp2(X,Y,lambda,x,y);
    cmatrix = zeros(n1,nr);
    cmatrix(1,:) =mu_lambda_interp;
    cmatrix(2,:) =0;
    cmatrix(3,:) =0;
    cmatrix(4,:) =mu_interp;
    cmatrix(5,:) =0;
    cmatrix(6,:) =lambda_interp;
    cmatrix(7,:) =mu_interp;
    cmatrix(8,:) =0;
    cmatrix(9,:) =0;
    cmatrix(10,:) =mu_interp;
    cmatrix(11,:) =lambda_interp;
    cmatrix(12,:) =0;
    cmatrix(13,:) =mu_interp;
    cmatrix(14,:) =0;
    cmatrix(15,:) =0;
    cmatrix(16,:) =mu_lambda_interp;        
end
end