function out = f(sigma,S0,K,r,q,T,market_price,flag)
% difference between BS model price and market price

d_1 = (log(S0/K)+(r - q +sigma^2/2)*T)/(sigma*sqrt(T));
d_2 = d_1 - sigma*sqrt(T);

if flag == 1   % put
    price_BS = -exp(-q*T)*S0*normcdf(-d_1)+K*exp(-r*T)*normcdf(-d_2);
else           % call
    price_BS = exp(-q*T)*S0*normcdf(d_1)-K*exp(-r*T)*normcdf(d_2);
end  

out = price_BS - market_price;
end

