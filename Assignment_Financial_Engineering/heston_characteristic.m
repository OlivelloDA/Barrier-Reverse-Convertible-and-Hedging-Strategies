function out = heston_characteristic(r, q, kappa, eta, theta, rho, sigma0, S0, u, t)

d = ((rho*theta*u*1i-kappa).^2-theta^2*(-1i*u-u.^2)).^0.5;
g = (kappa-rho*theta*u*1i-d)./(kappa-rho*theta*u*1i+d);

p1 = 1i.*u.*(log(S0)+(r-q).*t);
p2 =  eta.*kappa .*theta^(-2).*((kappa-rho*theta*u*1i-d).*t - ...
    2*log((1-g.*exp(-d*t))./(1-g)));
p3 = sigma0^2*theta^(-2).*(kappa-rho*theta*u*1i - d).*...
    (1- exp(-d*t))./(1-g.*exp(-d*t));

out = exp(p1).*exp(p2).*exp(p3);
end

