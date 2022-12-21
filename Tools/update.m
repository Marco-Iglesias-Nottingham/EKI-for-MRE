%% Marco Iglesias, Universtity of Nottingham, 2022
function Un=update(Un,Un_m,alpha,Z,Z_m,N_En)
%% Update of parameters via Kalman formula (in ensemble space)
% Note that this is different (but equivalent) to the update formula for the manuscript
% The formula from the paper is not efficient for large measurements

M=length(Z_m(:,1));
Delta_Z=sqrt(1/(N_En-1))*(Z-Z_m);

E=sqrt(alpha)*randn(M,N_En);  
E=E-mean(E,2);
Z=Z+E;

RHS=(Delta_Z')*Z;
B1=(1/alpha*(Delta_Z')*Delta_Z+eye(N_En))\RHS;
B2=Delta_Z'*(1/alpha*Z - 1/alpha^2*Delta_Z*B1);

for i=1:length(Un)
    Delta_U=sqrt(1/(N_En-1))*(Un{1,i}-Un_m{1,i});
    C_u_z=Delta_U;
    Un{1,i}=Un{1,i}-C_u_z*B2;
end




