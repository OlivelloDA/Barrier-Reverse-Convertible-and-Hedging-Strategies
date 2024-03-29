clear all; clc;

options = xlsread('Option_chain_15102021', 'options_cleaned');
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
sigma0 = 0.15;
kappa = 0.1;
eta = 0.1;
theta = 0.1; %higher than benchmark
rho = 0.4; %it should be between -1;1 so choose value in the middle

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
lb = [0.1,0.1,0.1,0.1,-1];
%fval is the value of the objective function, so in this case the error
%value
[error,fval]  =  fmincon(@(X)rmse_Heston(X,S0,K,r,q,T,market_price,flag),x0,A,B,Aeq,Beq,lb,ub);
%the only values that changes are those linked to X
disp(error)
disp(fval)

first = [   0.2020 ,  10.9242 ,   0.0840  ,  1.5515  ,  0.0135  ,  0.2191]
second = [0.3824  , 93.0843  ,  0.1252  , 59.6410 ,  -0.0092 , 0.4850]
third = [  0.2115    ,4.7020    ,0.1000    ,1.0439   , 0.0034  ,  0.2237]

sigma_opt = 0.2020 ;
kappa =  10.9242   ;
eta = 0.0840  ;
theta =  1.5515; %higher than benchmark
rho =   0.0135;


sigma_opt = 0.2195 ;
kappa =  0.3744   ;
eta = 0.4124  ;
theta =  0.6422; %higher than benchmark
rho =   0.0014;

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

% Monte Carlo properties : Reverse convertible at 6-months year Maturity

n = 129; % since the maturity is at 188 days --> 196 trading days , made through a proportion consideting 365 gg and 250 trading days
dt = 1/n;
rng('default')  %set the seed to let the simulations be consistent for delta 
m=50000;
T_exotic = 188/365
K_exotic = S0;
r_maturity = 0.000401829858818408; %already multiplied for the T_exotic
H = 0.80*S0;
K_reverse = S0 ;
S = zeros(m,n+1);
v = zeros(m,n+1);
S(:,1) = S0;
v(:,1) = sigma_opt^2;

% generate correlated random numbers
eps = normrnd(0,1,m,n);
epsS = normrnd(0,1,m,n);
eps1 = eps; 
eps2 = eps*rho + sqrt(1-rho^2)*epsS;

% simulate price paths according to Heston model

for j=2:n+1
        S(:,j) = S(:,j-1).*(1+(r_maturity-q)*dt+sqrt(v(:,j-1))*sqrt(dt).*eps1(:,j-1));
        v(:,j) = abs(v(:,j-1)+(kappa*(eta-v(:,j-1)))*dt+theta*sqrt(v(:,j-1))*...
            sqrt(dt).*eps2(:,j-1));  % reflection principle
end
figure()
plot(S)
%vanilla_call = exp(-r_maturity)*mean(max(S(:,end)-K_reverse,0))

  
% Price Down and out put barrier option

DOBP_dp = exp(-r_maturity).*max((H - min(S,[],2))./abs(H - min(S,[],2)), 0).*max(K_exotic-S(:,n+1),0);


%DOBP = mean(DOBP_dp)
%median(DOBP_dp)

%create indicator functions for BRC pricing
S_TH = zeros(m,1);
for j = 1:m
  
    if(S(j,n+1)<H)
        S_TH(j,1) = 1
    else 
    S_TH(j,1) = 0
    
    end
    
end

S_TS0 = zeros(m,1);

for j = 1:m
    for i = 1:n
  
        if(S(j,i)<=H)&(S(j,i+1)<S0)
        S_TS0(j,1) = S(j,n+1)/S0
        
    else 
    S_TS0(j,1) = 1
        end
    
    end
    
end

N = 1000
n_put = 10
c = 0.08
f = 0.045

price = (n_put*(K_exotic-S(:,n+1)).*S_TH + N.*S_TS0 +N*(c+f))-(N*exp(r_maturity)+ n_put*DOBP_dp)+N
t = exp(-r_maturity)*mean(price)



%HEDGING
%selling option
%using the money of the selling + borrowing money from the bank
%buy stocks for hedging

%we want to evaluate the option prices for every time step, the maturity
%changes at every time step. At the same time the stock changes in price,
%so the DELTA(t0,S0) value has changed in DELTA(t1,S1) and so I have to rebalance my portfolio.
%the cash flow coming from the rebalancing, is used to adjust the bank
%account. Assumes there are no "frctions" in rebalancing, no bid-ask
%spread, and buying/sell stuff at the same moment

%To calculate the Delta, use an approximation approach:
% Your second ‘approximation’ approach comes close.  First of all, you need to focus on the delta
% of the BRC. By definition, the delta of the BRC is the derivative 
% of the BRC’s price with respect to the price of the underlying (S0).
% That derivative can be approximated as
%     Delta » ( BRC(S0 + h) – BRC(S0 – h) ) / 2h
% which is related to your 2nd approach you compute the two prices in the numerator with 
% your Monte Carlo pricer. Make sure to obtain a stable approximation by fixing the seed of the 
% random numbers in the Monte Carlo simulation.
h = 0.55
S0 = 136
P1 = 999.0214
S0 = 134.95
P2 = 999.2228

Delta = (P1 - P2)/(2*h)