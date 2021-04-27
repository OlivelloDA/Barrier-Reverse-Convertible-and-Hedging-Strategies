function out = Heston_FFT(kappa, eta, theta, rho, sigma0, K, T, S0, r, q, flag)

%flag = 1 --> Call option
% flag = 0 --> Put option

% define parameters
N = 4096;         
alpha = 1.5;
eta_grid = 0.25;  
lambda = 2*pi/(N*eta_grid);
b = lambda*N/2;

% define grid of log-strikes
k = (-b:lambda:b-lambda); 

% compute rho
v = (0:eta_grid:(N-1)*eta_grid);
u = v-(alpha+1)*1i;
rho = exp(-r*T)* heston_characteristic(r,q,kappa,eta,theta,rho,sigma0,S0,u,T)./(alpha^2+alpha-v.^2+1i*(2*alpha+1)*v);


    simpson_1 = (1/3);                     
    simpson = ((3 + (-1).^(2:1:N))/3);
    simpson_int = [simpson_1 simpson];
    a = real(fft(rho.*exp(1i*v*b)*eta_grid.*simpson_int, N)); 

CallPrices = (1/pi)*exp(-alpha*k).*a;       

% find C(K,T)
KK = exp(k);
out = spline(KK,CallPrices,K); 

if flag == 0 
    out = out + K*exp(-r*T)- exp(-q*T)*S0;
end
end