clear all; clc;

options = xlsread('Option_chain_15102021', 'options');
stock = xlsread('Option_chain_15102021', 'Stock');
interest_rate = xlsread('Option_chain_15102021', 'interest_rate');

S0 = stock;
q = 0.02;
T = options(:,1);
K = options(:,2);
flag = options(:,3);
market_price = options(:,4);

% computation of the interest rate for each maturity 
r = spline(interest_rate(1,:), interest_rate(2,:)/100, T);
T = T/365; % maturity in years


% find the implied volatility
sigma = zeros(length(T),1);

for i=1:length(T)
    % f: difference 
    % @(sig): makes f a function with only one argument (sig)
   % sigma(i,1)= fzero(@(sig)f(sig,S0,K(i),r(i),q,T(i),market_price(i),flag(i)),[0,1]);
    % g: squared difference
    sigma(i,1)= fminbnd(@(sig)g(sig,S0,K(i),r(i),q,T(i),market_price(i),flag(i)),0,1);
%this functiom minimize tha value of g changing the sigma value.

end



%Price Vanilla Options under Heston Model using FFT and Carr-Madan formula

%Calibration Heston Characteristic function parameters

%define initial parameters(Rouah, 2013) 2*kappa*eta > theta^2(FellerCondition)
sigma0 = 0.7;
kappa = 0.2;
eta = 20;
theta = 0.30; %higher than benchmark
rho = -0.1; %it should be between -1;1 so choose value in the middle

%algorithm that optimizes the difference between the market price and the
%heston model price , giving upper and lower bound for every parameter
A = [];
B = [];
%A=[x(4)^2];
%b =[2*x(2)*x(3)];
x0 = [sigma0,kappa,eta,theta,rho];
Aeq = [];
Beq = [];
ub = [Inf,Inf,Inf,Inf,+1];
lb = [0,0,0,0,-1];

[error,fval]  =  fmincon(@(X)rmse_Heston(X,S0,K,r,q,T,market_price,flag),x0,A,B,Aeq,Beq,lb,ub);
%the only values that changes are those linked to X
disp(set)


second = [ 0.0553 ,   0.0181 , 167.6502   , 9.8462   , 0.9975,     0.8273]


% check the characteristic function
heston_characteristic(r,q,kappa,eta,theta,rho,sigma0,S0,1,1)
%last integer regards the integration rule: 0 = call , 1 = put(approximate
%carr-madan
Heston_FFT(kappa, eta, theta, rho, sigma0, K, T, S0, r, q, 0, 0) % rectangular rule
Heston_FFT(kappa, eta, theta, rho, sigma0, K, T, S0, r, q, 0, 1) % Simpson's rule


