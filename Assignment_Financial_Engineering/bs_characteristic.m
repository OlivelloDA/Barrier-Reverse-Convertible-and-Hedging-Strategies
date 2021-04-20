function out = bs_characteristic(r,q,sigma,S0,u,T)

out = exp(1i*u*(log(S0) + (r-q-sigma^2/2).*T)-(1/2)*sigma^2.*T.*u.^2);

end

