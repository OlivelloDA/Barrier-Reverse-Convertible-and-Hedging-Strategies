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
sigma0 = 0.3;
kappa = 0.2;
eta = 3;
theta = 0.7; %higher than benchmark
rho = 0.5; %it should be between -1;1 so choose value in the middle

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
%fval is the value of the objective function, so in this case the error
%value
[error,fval]  =  fmincon(@(X)rmse_Heston(X,S0,K,r,q,T,market_price,flag),x0,A,B,Aeq,Beq,lb,ub);
%the only values that changes are those linked to X
disp(error)
disp(fval)

first = [ 0.0605  ,  1.3936  , 10.0472  ,  4.7013 ,   0.5737 ,    1.8547e-05]
second = [0.2089   ,  0.1007   ,  1.0297  ,  0.4468  , -0.6557,   6.3826e-05]
third = [ 0.0356  ,   0.0738,    4.8792,   0.7911 ,  -0.0195,    3.0330e-05]
fourth = [  0.1355 ,   0.0837  ,  2.8938  ,  0.8324 ,  0.3967,    2.4504e-05]



sigma_opt = 0.0605 ;
kappa =  1.3936 ;
eta = 10.0472 ;
theta = 4.7013; %higher than benchmark
rho = 0.5737;

model_price = zeros(length(market_price),1)
for i=1:length(market_price)
    % f: difference 
    % @(sig): makes f a function with only one argument (sig)
   % sigma(i,1)= fzero(@(sig)f(sig,S0,K(i),r(i),q,T(i),market_price(i),flag(i)),[0,1]);
    % g: squared difference
    model_price(i,1)= Heston_FFT(kappa, eta, theta, rho, sigma_opt, K(i), T(i), S0, r(i), q, flag(i));
%this functiom minimize tha value of g changing the sigma value.

end

figure()
plot(K,model_price,'r*','LineWidth',1.1)
hold on
plot(K,market_price,'bo','LineWidth',1.1)
hold on
xlabel('K')
ylabel('Price')
title('Calibration: Heston')
legend('Heston Price','Market Price','Market')
disp(['RMSE: ' num2str(sqrt(sum((model_price-market_price)).^2)/length(market_price))]);

% Monte Carlo properties
m=100000;        % number of price paths
n=T_exotic_day;  % number of time steps: 47
dt = 1/365;      % length of time step (daily)


S = zeros(m,n+1);
Z = normrnd(0,1,m,n); % sample random numbers
S(:,1) = S0;

% fill the matrix with Formula (1)
for j=2:n+1
    S(:,j) = S(:,j-1).*exp((r_exotic-q-0.5*sigma^2)*dt + sigma*sqrt(dt)*Z(:,j-1));
end


% a) Asian call option
mean_S = mean(S,2); % mean of the stock price along each sample path
AC_dp = exp(-r_exotic*T_exotic).*max(mean_S-K_exotic, 0);
AC = mean(AC_dp)

% b) barrier options
UIBP_dp = exp(-r_exotic*T_exotic).*max((max(S,[],2)-H)./abs(max(S,[],2)-H), 0).*max(K_exotic-S(:,n+1),0);
UIBP = mean(UIBP_dp)

UOBP_dp = exp(-r_exotic*T_exotic).*max((H - max(S,[],2))./abs(H - max(S,[],2)), 0).*max(K_exotic-S(:,n+1),0);
UOBP = mean(UOBP_dp)

