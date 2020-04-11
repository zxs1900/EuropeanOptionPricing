clc;
clear all;
close all;

%%
%Parameters

S0 = 100;    %initial stock price
K = 110;     %strike price
r = 0.05;    %interest rate 
T = 3/12;    %time to maturity
sigma = 0.2; %volatility
x_min = log(S0) - 2 * sigma * sqrt(T); %maximum stock price
x_max = log(S0) + 2 * sigma * sqrt(T); %minimum stock price


N = 100; %time steps
dt = T/N;
M = 30; %price steps
dx = (x_max - x_min)/M;
N = round(T/dt);
dt = T/N;



EUput = zeros(M+1,N+1); %grid setup
Amput = zeros(M+1,N+1);
EUEXput = zeros(M+1,N+1);
AMEXput = zeros(M+1,N+1);
vetx = linspace(x_min,x_max,M+1)';
vetS = exp(vetx);
vetj = 0:N;

% set up coefficients 
a  = 0.5*dt*((sigma^2/dx^2) - (r - 0.5 * sigma^2) / dx);
b = 1 - dt * sigma^2 / dx^2;
c  =0.5*dt*((sigma^2/dx^2) + (r - 0.5 * sigma^2) / dx);

%%
%European put option

% set up boundary conditions
EUput(:,N+1) = max(K-vetS,0);
EUput(1,:)=EUput(1,N+1)*exp(-r*dt*(N-vetj));
EUput(M+1,:) = 0;

% solve backward in time
for j=N:-1:1
   for i=2:M
      EUput(i,j) = (a*EUput(i-1,j+1) + b*EUput(i,j+1)+ ...
         c*EUput(i+1,j+1)) * (1/(1 + r*dt));
   end
end
% return price, possibly by linear interpolation outside the grid
p_EUput = interp1(vetS, EUput(:,1), S0);

%%
%American put option

% set up boundary conditions
AMput(:,N+1) = max(K-exp(vetx),0);
AMput(1,:) = AMput(1,N+1);
AMput(M+1,:) = 0;

ab = zeros(M+1,1);
ab = max(K-exp(vetx),0);
ab = ab(2:M);

% solve backward in time
for j=N:-1:1
   for i=2:M
       
      AMput(i,j) = (a*AMput(i-1,j+1) + b*AMput(i,j+1)+ ...
         c*AMput(i+1,j+1)) * (1/(1 + r*dt));
     AMput(i,j) = max(AMput(i,j),ab(i-1));
   end
end
% return price, possibly by linear interpolation outside the grid
p_AMput = interp1(vetS, AMput(:,1), S0);


%figure(2);
plot(vetS,EUput(:,1),vetS, AMput(:,1),vetS, AMput(:,end))
legend({'European put','American put','Payoff'});
title('Vanilla options');

%%
%Eu Exotic put option

% set up boundary conditions
EUEXput(:,N+1) = max(K - 100 *(vetS/S0).^(1/T),0);
EUEXput(1,:) =EUEXput(1,N+1)*exp(-r*dt*(N-vetj));
EUEXput(M+1,:) = 0;

% solve backward in time
for j=N:-1:1
   for i=2:M
      EUEXput(i,j) = (a*EUEXput(i-1,j+1) + b*EUEXput(i,j+1)+ ...
         c*EUEXput(i+1,j+1)) * (1/(1 + r*dt));
   end
end
% return price, possibly by linear interpolation outside the grid
p_EUEXput = interp1(vetS, EUEXput(:,1), S0);





%%
%Am Exotic put option


% set up boundary conditions

AMEXput(:,N+1) = max(K - S0 *(vetS/S0).^(1/T),0);
AMEXput(1,:) = AMEXput(1,N+1);
AMEXput(M+1,:) = 0;
% set up coefficients 
a  = 0.5*dt*((sigma^2/dx^2) - (r - 0.5 * sigma^2) / dx);
b = 1 - dt * sigma^2 / dx^2;
c  =0.5*dt*((sigma^2/dx^2) + (r - 0.5 * sigma^2) / dx);

ab = zeros(M+1,1);
ab = max(K-exp(vetx),0);
ab = ab(2:M);

% solve backward in time
for j=N:-1:1
   for i=2:M
       
      AMEXput(i,j) = (a*AMEXput(i-1,j+1) + b*AMEXput(i,j+1)+ ...
         c*AMEXput(i+1,j+1)) * (1/(1 + r*dt));
     AMEXput(i,j) = max(AMEXput(i,j),ab(i-1));
   end
end
% return price, possibly by linear interpolation outside the grid
p_AMEXput = interp1(vetS, AMEXput(:,1), S0);

figure();
plot(vetS, EUEXput(:,1),vetS, AMEXput(:,1),vetS, AMEXput(:,end))
legend({'European put','American put','Payoff'});
title('Exotic options');

