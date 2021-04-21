function out = rmse_Heston(X,S0,K,r,q,T,market_price,flag)
% computes the sum of squared differences (errors) for a vector of options
model_price = zeros(length(market_price),1);

sigma0 = X(1);

kappa = X(2);

eta = X(3);

theta = X(4);

rho = X(5);

for i=1:length(market_price)
    % f: difference 
    % @(sig): makes f a function with only one argument (sig)
   % sigma(i,1)= fzero(@(sig)f(sig,S0,K(i),r(i),q,T(i),market_price(i),flag(i)),[0,1]);
    % g: squared difference
    model_price(i,1)= Heston_FFT(kappa, eta, theta, rho, sigma0, K(i), T(i), S0, r(i), q, flag(i));
%this functiom minimize tha value of g changing the sigma value.

end

% RMSE for a set of parameters applied to all the options market prices
out = sqrt(sum( model_price - market_price)/length(market_price));

end

