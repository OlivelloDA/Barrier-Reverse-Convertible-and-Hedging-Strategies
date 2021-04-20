function out = BS_FFT(sigma, K, T, S0, r, q, type, integration_rule)
% integration_rule = 0 --> rectangular rule
% integration_rule = 1 --> Simpson's rule
% type = 0 --> call option
% type = 1 --> put option

% define parameters
N = 4096;         
alpha = 1.5;
eta = 0.25;  
lambda = 2*pi/(N*eta);
b = lambda*N/2;

% define grid of log-strikes
k = (-b:lambda:b-lambda); 

% compute rho
v = (0:eta:(N-1)*eta);
u = v-(alpha+1)*1i;
rho = exp(-r.*T).*bs_characteristic(r, q, sigma, S0, u , T)./(alpha^2+alpha-v.^2+1i*(2*alpha+1)*v);

if integration_rule == 0                    % rectangular rule
    a = real(fft(rho.*exp(1i*v*b)*eta, N)); 

elseif integration_rule == 1                % Simpson's rule
    simpson_1 = (1/3);                     
    simpson = ((3 + (-1).^(2:1:N))/3);
    simpson_int = [simpson_1 simpson];
    a = real(fft(rho.*exp(1i*v*b)*eta.*simpson_int, N)); 
end
CallPrices = (1/pi)*exp(-alpha*k).*a;       

% find C(K,T)
KK = exp(k);
out = spline(KK,CallPrices,K); 

% if put option
if type ==1 
    out = out + K*exp(-r*T)- exp(-q*T)*S0;
end
end