% Exercise 1
clear all; clc;

S0 = 100;
K = 105; 
r = 0.04; 
T = 2;
sigma = 0.25;
q = 0.05;

% Price the call option
call =  BS_call(sigma, S0, K, r, q, T)

% Price the put option
put =  BS_put(sigma, S0, K, r, q, T)

% Price the put option (with the put-call parity)
put_pcp = -exp(-q*T)*S0 + K*exp(-r*T) + call