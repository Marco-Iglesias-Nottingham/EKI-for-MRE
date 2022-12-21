%% Marco Iglesias, Universtity of Nottingham, 2022
clear 
close all
addpath('Tools')
%% This code runs the 2D implementation of ensemble Kalman Inversion (EKI)
% that we used for some of the 2D examples in Marco Iglesias et al 2022 Phys. Med. Biol. 67 235003




%% Generate Synthetic data

rng(426922); %%seed for random numbers

%call the Truth and compute noise-free displacements on a M x M grid
load Truth;  
data_clean=Solver(Truth,reshape(Truth.Storage,Truth.Grid.Nx,Truth.Grid.Ny),...'
    reshape(Truth.Loss,Truth.Grid.Nx,Truth.Grid.Ny));
%call noisy data
sigma = 0.025;
[data_noisy, inv_sqrt_Gamma]=get_Data(data_clean,sigma,0.0005);
%plot noisy data
plot_data(data_noisy, Truth)


%% Define forward model for inversion
Nx = 120;	% NxN square grid foir the inversion
hmax=0.0009;  % mesh size for PDE tool
M=length(Truth.Grid.X_meas);% get number of measurements (grid M xM)
Model=set_model(Nx,hmax,M); % here we hard coded model parameters


%% Define prior
prior.no_fields=4; % random fields needed 
% for stifness properties in medulla, cortex, background and cyst.

% hyper parameters for Level-set function
level.nu = 3.0; level.sigma = 1.0;
level.len(1).mean.lim=[0.03,0.1]; 
level.len(2)=level.len(1);

% hyper parameters for Storage and loss
% Storage and loss log-means (uniform priors on the mean and here we 
% specify end points of interval)
Storage(1).mean.lim=[log(3700),log(4600)];
Storage(2).mean.lim=[log(2450),log(3800)];
Storage(3).mean.lim=[log(900),log(1700)];
Storage(4).mean.lim=[log(150),log(400)];
Loss(1).mean.lim=[log(2000),log(2900)];
Loss(2).mean.lim=[log(2200),log(2900)];
Loss(3).mean.lim=[log(1100),log(2000)];
Loss(4).mean.lim=[log(3100),log(3500)];

% signal variance (fixed), smoothness parameter (fixed) and 
% lengthscales (uniform priors for these and here we specify end points of
% interval
Storage(1).sigma = 0.01;
Storage(2).sigma = 0.03;
Storage(3).sigma = 0.03;
Storage(4).sigma = 0.1;
for i=1:prior.no_fields
    Storage(i).nu = 1.0; 
    Loss(i).nu = 1.0; 
    Storage(i).len(1).lim=[0.035,0.3];
    Loss(i).len(1).lim=[0.035,0.3];
    Storage(i).len(2).lim=[0.035,0.3];
    Loss(i).len(2).lim=[0.035,0.3];
    Loss(i).sigma=Storage(i).sigma;
end
    
Loss(4).sigma = 0.01;
prior.level=level; prior.Storage=Storage; prior.Loss=Loss;

% we work on a transform for the lengthscales to ensure these are always
% in the interval [0.015, 0.4] (this is relative to the unit interval).
prior.max_len=0.4;
prior.min_len=0.015;

% threshold for level-set function (we tried this random but did not have
% any massive effect
prior.threshold=1.75;

%% Inversion

N_En=500; %Ensemble size  
% (Note we use larger ensemble in some examples of the paper and we ran
% those in HPC)

rng(989*2)

% Sample from the prior in the parameter space
Un=Get_prior_fields(prior,Model,N_En);

load mask; %get mask
Model.mask=mask;
% run Ensemble Kalman Inversion (EKI)
EKI(Un,Model,prior,N_En,data_noisy,inv_sqrt_Gamma);

% plot prior/posterior cyst probability and mean/variances for Storage/Loss 
% modulus 
plot_mean_variance;







