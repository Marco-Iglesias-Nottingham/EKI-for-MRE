%% Marco Iglesias, Universtity of Nottingham, 2022
function [data, inv_sqrt_Gamma]=get_Data(data_clean,sigma1,sigma2)
%% add noise to the noise-free displacements as described in:
% Marco Iglesias et al 2022 Phys. Med. Biol. 67 235003
% this also returns the inverse of the square root of the measurement error
% covariance (which we assume is diagonal)
N=length(data_clean)/4;
sqrt_Gamma=zeros(N*4,1);
data=zeros(N*4,1);
for i=1:4
    d1=data_clean(1+(i-1)*N:i*N);
    M1=length(d1);
    Gamma1_d1=(sigma1*abs(d1)).^2;
    Gamma2_d1=(sigma2*abs(max(d1)-min(d1))).^2;
    Gamma_d1=Gamma1_d1+Gamma2_d1;
    sqrt_Gamma(1+(i-1)*N:i*N)=sqrt(Gamma_d1);
    noise_d1 =sqrt(Gamma1_d1).*normrnd(0,1,M1,1)+sqrt(Gamma2_d1).*normrnd(0,1,M1,1);
    data(1+(i-1)*N:i*N)=d1+noise_d1;
end
inv_sqrt_Gamma=1./sqrt_Gamma;

