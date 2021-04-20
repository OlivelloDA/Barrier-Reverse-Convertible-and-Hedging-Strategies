function out = rmse_Heston(sigma0,kappa,S0,K,r,q,T,market_price,flag,eta,theta,rho)
% computes the sum of squared differences (errors) for a vector of options
model_price = zeros(length(market_price),1);
for i=1:length(market_price)
    % f: difference 
    % @(sig): makes f a function with only one argument (sig)
   % sigma(i,1)= fzero(@(sig)f(sig,S0,K(i),r(i),q,T(i),market_price(i),flag(i)),[0,1]);
    % g: squared difference
    model_price(i,1)= Heston_FFT(kappa, eta, theta, rho, sigma0, K(i), T(i), S0, r(i), q, flag(i), flag(i));
%this functiom minimize tha value of g changing the sigma value.

end

% vector of Heston prices
out = sqrt(sum( model_price - market_price)/length(market_price));

end


