%% Exercise 3

K_exotic=1500;
H=1570;
T_exotic_day=47;
T_exotic=T_exotic_day/365;

interest_rate = xlsread('SP500_03_01_00', 'interest_rate');
stock = xlsread('SP500_03_01_00', 'stock');
S0 = stock;
q=0;

% computation of the interest rate
r_exotic = spline(interest_rate(:,1), interest_rate(:,2)/100, T_exotic_day);

% use the volatility (sigma) from calibration
sigma=0.2228;

% Monte Carlo properties
m=100000;        % number of price paths
n=T_exotic_day;  % number of time steps: 47
dt = 1/365;      % length of time step (daily)


S = zeros(m,n+1);
Z = normrnd(0,1,m,n); % sample random numbers
S(:,1) = S0;

% fill the matrix with Formula (1)
for j=2:n+1
    S(:,j) = S(:,j-1).*exp((r_exotic-q-0.5*sigma^2)*dt + sigma*sqrt(dt)*Z(:,j-1));
end


% a) Asian call option
mean_S = mean(S,2); % mean of the stock price along each sample path
AC_dp = exp(-r_exotic*T_exotic).*max(mean_S-K_exotic, 0);
AC = mean(AC_dp)

% b) barrier options
UIBP_dp = exp(-r_exotic*T_exotic).*max((max(S,[],2)-H)./abs(max(S,[],2)-H), 0).*max(K_exotic-S(:,n+1),0);
UIBP = mean(UIBP_dp)

UOBP_dp = exp(-r_exotic*T_exotic).*max((H - max(S,[],2))./abs(H - max(S,[],2)), 0).*max(K_exotic-S(:,n+1),0);
UOBP = mean(UOBP_dp)
