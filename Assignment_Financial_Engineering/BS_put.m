function out = BS_put(sigma, S0, K, r, q, T) 
% Calculates the price of a European put option

d_1 = (log(S0/K)+(r - q +sigma^2/2)*T)/(sigma*sqrt(T));
d_2 = d_1 - sigma*sqrt(T);

out = K*exp(-r*T)*normcdf(-d_2)-S0*exp(-q*T)*normcdf(-d_1);
end 