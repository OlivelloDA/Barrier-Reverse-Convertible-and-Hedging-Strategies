%% Exercise session 3 

r = 0.05;
q = 0.01;
T = 1;
K = [90,100,110];
S0 = 100;
kappa = 0.5;
theta = 0.2;
eta = 0.05;
v0 = 0.05;
rho = -0.75;

%% Exercise 1 
method = 1; 
% method = 1 -> Euler
% method = 0 -> Milstein

n = 252*T; % daily steps
dt = T/n;  

m=100000;

S = zeros(m,n+1);
v = zeros(m,n+1);
S(:,1) = S0;
v(:,1) = v0;

% generate correlated random numbers
eps = normrnd(0,1,m,n);
epsS = normrnd(0,1,m,n);
eps1 = eps; 
eps2 = eps*rho + sqrt(1-rho^2)*epsS;

% simulate price paths
if method == 1  % Euler
    for j=2:n+1
        S(:,j) = S(:,j-1).*(1+(r-q)*dt+sqrt(v(:,j-1))*sqrt(dt).*eps1(:,j-1));
        v(:,j) = abs(v(:,j-1)+(kappa*(eta-v(:,j-1)))*dt+theta*sqrt(v(:,j-1))*...
            sqrt(dt).*eps2(:,j-1));  % reflection principle
    end
else  % Milstein 
    for j=2:n+1
        S(:,j) = S(:,j-1).*(1+(r-q)*dt+sqrt(v(:,j-1))*sqrt(dt).*eps1(:,j-1));
        v(:,j) = abs(v(:,j-1)+ (kappa*(eta-v(:,j-1))-theta^2/4)*dt+theta*sqrt(v(:,j-1))*sqrt(dt).*eps2(:,j-1)+theta^2*dt.*(eps2(:,j-1).^2)/4);
    end
end

MC = exp(-r*T)*mean(max(S(:,end)-K,0))

% compute the price with FFT
FFT = Heston_FFT(kappa, eta, theta, rho, sqrt(v0), K, T, S0, r, q, 0,1)


%% Exercise 2

H = 90;

% digital down-and-out option
MC = exp(-r*T)*mean(max(0, (min(S,[],2)-H)./abs(min(S,[],2)-H)))


