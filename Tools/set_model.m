%% Marco Iglesias, Universtity of Nottingham, 2022
function Model=set_model(N,hmax,M)

%%%%%%%%%%%%%%%%%%%%%%
% Define parameters for the 2D forward model and random fields. This is all
% saved in structure Model

%% Here we define the grid for random fields of size Nx x Nx
Nx=N;
Ny=Nx;
Grid.dim=2;
Grid.Nx=Nx; 
Grid.Ny=Ny; 
Grid.N=Grid.Nx*Grid.Ny;
Model.Grid=Grid;

%% Parameters for viscoelastic model

Model.rho=1000; %density
Model.omega=60*2*pi; %angular velocity
Model.dispAmplitude=1e-6;  %displacement amplitude on boundary, 
Model.nu=0.499; %Poisson's ratio
Lx=0.12;       %edges of the square (in meters)
Ly=Lx;

[X,Y]=meshgrid(linspace(0,Lx,Nx),linspace(0,Ly,Nx));
Model.Grid.X=X;
Model.Grid.Y=Y;

% Grid for measurements M x M
[Model.Grid.X_meas,Model.Grid.Y_meas]=meshgrid(linspace(0,Lx,M),linspace(0,Ly,M));


%% Below is stuff needed for the PDE tool box (e.g. mesh generation) so that we do not 
% define this everytime we run the forward model
lowerLeft  = [0   ,0   ];
lowerRight = [Lx , 0  ];
upperRight = [Lx , Ly];
upperLeft =  [0.0 , Ly];
% Geometry matrix
S = [3,4 lowerLeft(1), lowerRight(1), upperRight(1), upperLeft(1), ...
    lowerLeft(2), lowerRight(2), upperRight(2), upperLeft(2)];
gdm = S';
% Names
ns = 'S';
% Set formula
sf = 'S';
% Invoke decsg
g = decsg(gdm,ns,sf');
% Import g into model using geometryFromEdges.
model = createpde(2);
geometryFromEdges(model,g);
generateMesh(model,'Hmax',hmax,'GeometricOrder','linear');
%generateMesh(model,'Hmax',0.001,'GeometricOrder','linear');

%pdemesh(model)
Model.Lx=Lx;
Model.Ly=Ly;
Model.FM=model;


