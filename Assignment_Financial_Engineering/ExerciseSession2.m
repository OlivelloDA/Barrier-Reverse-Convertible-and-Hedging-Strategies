%% Exercise 1
S0 = 100;
q=0;
r=0.05; 
sigma= 0.2;
u=1;
T=1;
bs_characteristic(r,q,sigma,S0,u,T)



%% Exercise 2
clear;

S0 = 100; 
q = 0.01; 
r = 0.05;
sigma = 0.25;
K = [90,95,100,105,110];
T = 3/12;

BS_call(sigma, S0, K, r, q, T)       % closed form 
BS_FFT(sigma, K, T, S0, r, q, 0, 0)  % FFT with rectangular rule
BS_FFT(sigma, K, T, S0, r, q, 0, 1)  % FFT with Simpson's rule


%% Exercise 3
clear;

S0 = 100; 
q = 0.01; 
r = 0.05;
K = [90,95,100,105,110];
T = 3/12;

v0 = 0.05; sigma0 = sqrt(v0);
kappa = 0.5;
eta = 0.05;
theta = 0.2;
rho = -0.75;

% check the characteristic function
heston_characteristic(r,q,kappa,eta,theta,rho,sigma0,S0,1,1)

Heston_FFT(kappa, eta, theta, rho, sigma0, K, T, S0, r, q, 0, 0) % rectangular rule
Heston_FFT(kappa, eta, theta, rho, sigma0, K, T, S0, r, q, 0, 1) % Simpson's rule
