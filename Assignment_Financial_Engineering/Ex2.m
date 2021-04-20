%% Exercise 2a
clear all; clc;

option = xlsread('SP500_03_01_00', 'option');
stock = xlsread('SP500_03_01_00', 'stock');
interest_rate = xlsread('SP500_03_01_00', 'interest_rate');

S0 = stock;
q = 0;
T = option(:,1);
K = option(:,2);
flag = option(:,3);
market_price = option(:,4);

% computation of the interest rate for each maturity 
r = spline(interest_rate(:,1), interest_rate(:,2)/100, T);
T = T/365; % maturity in years


% find the implied volatility
sigma = zeros(length(T),2);

for i=1:length(T)
    % f: difference 
    % @(sig): makes f a function with only one argument (sig)
    %sigma(i,1)= fzero(@(sig)f(sig,S0,K(i),r(i),q,T(i),market_price(i),flag(i)),[0,1]);
    % g: squared difference
    sigma(i,2)= fminbnd(@(sig)g(sig,S0,K(i),r(i),q,T(i),market_price(i),flag(i)),0,1);
%this functiom minimize tha value of g changing the sigma value.
end

disp(sigma(1:6,:)); 


%% Exercise 2b

TT=unique(T);

figure()
for i=1:length(TT)
   KK=K(T==TT(i));
   vol=sigma(T==TT(i))';
   plot(KK,vol, 'LineWidth',i)
   hold on
end

xlabel('K')
ylabel('\sigma(K,T)')
title('Implied Volatility')
legend(strcat('T=',num2str(TT,3)))

%% Exercise 2c

sigma_opt = fminbnd(@(sig)sse(sig,S0,K,r,q,T,market_price,flag),0,1);

% plot the model prices and market prices against the strikes
figure()
plot(K,BS_price(sigma_opt,S0,K,r,q,T,flag),'r*','LineWidth',1.1)
hold on
plot(K,BS_price(sigma(:,1),S0,K,r,q,T,flag),'g+','LineWidth',1.1)
hold on
plot(K,market_price, 'o','LineWidth',1.1)
hold off
xlabel('K')
ylabel('Price')
title('Calibration: Black-Scholes')
legend('BS(\sigma)','BS(\sigma(K,T))','Market')

