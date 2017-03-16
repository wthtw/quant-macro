% This program solves the stochastic growth model by value function
% iteration. "simulate_stochastic_growth.m" should be in the same
% directory as the "stochastic_growth_model.m". You should run
% "stochastic_growth_model.m"

clf
clear
NK=101;

sigmaR = 1; % Relative risk aversion coefficient

beta  = 0.987;
alpha = 0.4;
delta = 0.012;

rho   = 0.95;
sigma = 0.007;

NZ   = 5; % number of "z" states
mz   = 3; % the bounds of the grid is chosen to be -m*sigma and m*sigma
zz(1) = -mz*sigma/(sqrt(1-rho^2));
zz(NZ)=  mz*sigma/(sqrt(1-rho^2));
d     = (zz(NZ)-zz(1))/(NZ-1); % distance between grids

% form the grid space for the shock "z"
for i=2:NZ-1
    zz(i) = zz(i-1) + d;
end
z = exp(zz);

% form the transition probabilities
for i=1:NZ
    pi(i,1)  =     normcdf((zz(1) +d/2-rho*zz(i))/sigma,0,1);    
    pi(i,NZ) = 1 - normcdf((zz(NZ)-d/2-rho*zz(i))/sigma,0,1);
    for j=2:(NZ-1)
        pi(i,j) = normcdf((zz(j)+d/2-rho*zz(i))/sigma,0,1)  ...
                - normcdf((zz(j)-d/2-rho*zz(i))/sigma,0,1);
    end
end

% compute the cumulative distribution function from the transition matrix
for iz=1:NZ
    for jz=1:NZ 
        cdfz(iz,jz) = sum(pi(iz,1:jz));
    end 
end

  
% Compute the steady state capital stock and form the grid space
kss_min = ((1/beta - 1 + delta)/(alpha*z(1)))^(1/(alpha-1));
kss_max = ((1/beta - 1 + delta)/(alpha*z(NZ)))^(1/(alpha-1));
kl = kss_min;
kh = kss_max;
K = [kl : (kh - kl )/(NK-1) : kh];

v(1:NZ,1:NK) = 0;
distance = 1;
while distance > 0.00001
    % solve tha value for all z and k
    for i=1:NZ
    for j=1:NK
        
        U = log(((1-delta)*K(j)+z(i)*K(j).^alpha-K).*((1-delta)*K(j)+z(i)*K(j).^alpha-K>0))...
            + beta * pi(i,:) * v;

%          if ( (sigmaR > 1) | (sigmaR<1) )
%           U = (((1-delta)*K(j)+z(i)*K(j).^alpha-K).*((1-delta)*K(j)+z(i)*K(j).^alpha-K>0) +eps).^(1-sigmaR)/(1-sigmaR)...
%             + beta * pi(i,:) * v;
%          end

        [tv(i,j)  indxg(i,j) ] = max(U); 
   
    end
    end
    distance = max(max(abs(v-tv)));
    v=tv;         
end


% compare the the decision rule with the ones
% corresponding to analytical solution for delta=1
clf
if delta==1    
    for i=1:NZ
        g_a(i,:) = alpha * beta * z(i) * K.^alpha;        
        % draw the analytical solution
        plot(K, g_a(i,:),'-')    
        hold on
        plot(K, K(indxg(i,:)),'.')
        %pause
    end
    % plot the 45 degree line
    plot(K,K)
    pause    
end

% run the program that simulates the economy given the decision rules above
% note that the file "simulate_stochastic_growth.m" should be in the same
% directory as the "stochastic_growth_model.m"

Simulate_Stochastic_growth_model