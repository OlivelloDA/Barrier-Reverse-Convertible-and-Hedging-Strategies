function out = BS_price(sigma,S0,K,r,q,T,flag) 
% faster method to compute Black-Scholes prices of call and put options

d_1 = (log(S0./K)+(r - q +sigma.^2/2).*T)./(sigma.*sqrt(T));
d_2 = d_1 - sigma.*sqrt(T);

call = exp(-q.*T).*S0.*normcdf(d_1)-K.*exp(-r.*T).*normcdf(d_2);
put = -exp(-q.*T).*S0.*normcdf(-d_1)+K.*exp(-r.*T).*normcdf(-d_2);

% put option: flag = 1
% call option: flag = 0
out = call.*(1-flag) + put.*flag;
end

