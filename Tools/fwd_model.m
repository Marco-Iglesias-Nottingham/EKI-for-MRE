%% Marco Iglesias, Universtity of Nottingham, 2022
function [out, Storage,Loss,Ind]=fwd_model(Un,Model,prior)
%% This function computes the parameter-to-output map. 
% It takes the parameters in Un and 
% computes the physical properties (Storage and Loss) 
Grid=Model.Grid;
[Storage, Loss,Level]=physical(Un,Model,prior);
% which are then needed to run the viscoleastic model
% and produce predictions on measurement grid:
out=Solver(Model,reshape(Storage,Model.Grid.Nx,Grid.Ny),reshape(Loss,Model.Grid.Nx,Grid.Ny));
% We also compute the indicator function of the cyst to then visualise
% probability of cyst
Ind=zeros(Model.Grid.N,1);
Ind(Level>=prior.threshold)=1;
