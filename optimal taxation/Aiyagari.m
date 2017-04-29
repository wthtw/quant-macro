clear

% parameters
%%%%%%%%%%%%
z = 1;
alpha = 0.4;
delta = 0.08;
beta  = 0.96;
rho   = 0.90; % parameters of labor endowment shock process l(t) = "phi(t)"
sigma = 0.20; % log(l(t+1)) = rho * log(l(t)) + epsilon(t) ; 
              % epsilon(t) is normally distributed with mean zero and
              % standard deviation sigma.
tau_y_bench = 0.3;
tau_c_bench = 0.075;
lambda = 2;
mu = 0.10;

NA = 400;   % number of intervals in A grid-space (analogous to K).
NL  = 5; % number of "l" states (analgous to Z).
%%%%%%%%%%%%

% Form the grid space
al = 0;
ah = 25;
A = [al : (ah - al )/(NA-1) : ah];


% Here we find a Markov chain which approximates the AR(1) process
% log(l(t+1)) = rho * log(l(t)) + epsilon(t), with mean 0 and s.d. sigma.
% We store this in the transition matrix "pil".

m   = 3; % the bounds of the grid is chosen to be -m*sigma and m*sigma
ll(1) = -m*sigma/(sqrt(1-rho^2));
ll(NL)=  m*sigma/(sqrt(1-rho^2));
d     = (ll(NL)-ll(1))/(NL-1); % distance between grids

% form the grid space for the shock "l" 
for i=2:NL-1
    ll(i) = ll(i-1) + d;
end
l = exp(ll);

% form the transition probabilities
for i=1:NL
    pil(i,1)  =     normcdf((ll(1) +d/2-rho*ll(i))/sigma,0,1);    
    pil(i,NL) = 1 - normcdf((ll(NL)-d/2-rho*ll(i))/sigma,0,1);
    for j=2:(NL-1)
        pil(i,j) = normcdf((ll(j)+d/2-rho*ll(i))/sigma,0,1)  ...
                 - normcdf((ll(j)-d/2-rho*ll(i))/sigma,0,1);
    end
end
pil

% here we form the stationary distribution and the value function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PIL  = [1 0 0 0 0];
PIL  = PIL* pil^1000; % stationary distribution of agents over "l"
LBAR = PIL * l';      % mean labor endowment

% make the initial guess for the value function
V_benchmark(1:NL,1:NA) = 0;

% make an initial guess for the stationary distribution
psi(1:NL,1:NA) = 0;
psi(1,1)       = 1;

% make an initial guess for interest rate
dist_r = 1;
r = 1/beta - 1 - 0.001;
r = 0.0379;

while dist_r > 0.0001
    R = (1 + (r * (1 - tau_y_bench)));
    KBAR = LBAR * ((r+delta)/(z*alpha)) ^(1/(alpha-1));
    w = (1-alpha)*(KBAR/LBAR)^alpha;
  
    % Solve the consumer's problem with value function iteration
    dist_V = 1; % this is our starting threshhold for the error
    while dist_V > 0.0000001    % we specify a convergence criterion
        for i=1:NL
            for j=1:NA
                %%%%%% THE CODE FOR N AND NSTAR ISN'T CORRECT%%%%%% %%%%%% 
                n = (1/(1 + lambda)) - ((1/(1 + lambda))*((R*A(j) - A)/((1 - tau_y_bench)*w*l(i))));
                n_star = min(max(n, 0), 1);
                %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% 
                %C = (1+r)*A(j) + w * l(i) - A;
                C = (R /(1 + tau_c_bench))*A(j) + (((1 - tau_y_bench)/(1 + tau_c_bench))* w * l(i) * n_star) - ((1/(1 + tau_c_bench)) * A);    % EDITED
                U = log(C.*(C>0)) + (lambda * log(1 - n_star)) + beta * pil(i,:) * V_benchmark;   % EDITED
                % note: do not replace n_star above with n, you will get
                % NaN errors.
                [TV(i,j) indxg(i,j)] = max(U);
            end
        end
        dist_V = max(max(abs(V_benchmark-TV)));
        V_benchmark=TV;
    end
    %%%%%%%%%% Consumer's problem has been solved.
    
    % here we find the stationary distribution of agent 
    dist_psi = 1;   % this is our starting threshold for the error
    while dist_psi > 0.0000000001   % we specify a convergence criterion
        for ip=1:NL   % lprime
            for jp=1:NA    % aprime 

               temp_sum=0;
               for i=1:NL   % summing over l
                   for j=1:NA   % summing over a
                       temp_sum = temp_sum + psi(i,j) * pil(i,ip) * (indxg(i,j) == jp);
                   end
               end
               
               Tpsi(ip,jp) = temp_sum;
            end
        end
        dist_psi = max(max(abs(psi-Tpsi))) ;
        psi=Tpsi;
    end
    % psi found.
    
    KBARn = 0;
    LBARn = 0;
    for i=1:NL  % l
        for j=1:NA  % a
            KBARn = KBARn + A(j) * psi(i,j);
            LBARn = LBARn + l(i) * psi(i,j);
        end
    end
    % computed K' and L'
    
    rn = z*alpha*(KBARn/LBAR)^(alpha-1) -delta;
    [r rn]
    dist_r = abs(r-rn)
    
    % update interest rate
    r = (1 - mu) * r + mu * rn;
    
    % computed r'
end


optimal_n = n_star(indxg);
optimal_c = A(indxg);

% K_bench = 0;
% for i = 1:NL
%     for j = 1:NA
%         K_bench = K_bench + A(j)*psi(i,j);
%     end
% end
K_bench = KBARn;

% N_bench = 0;
% for i = 1:NL
%     for j = 1:NA
%         N_bench = N_bench + l(i)*psi(i,j);
%     end
% end
N_bench = LBARn;

Y_bench = K_bench.^(alpha)*(N_bench.^(1 - alpha));

aggregate_c_bench = 0;
for i = 1:NL
    for j = 1:NA
        aggregate_c_bench = aggregate_c_bench + optimal_c(i,j)*psi(i,j);
    end
end
% notice that aggregate_c == K, meaning that we have equilibrium in the
% goods market.

G_bench = (tau_y_bench * (Y_bench - delta * K_bench)) + (tau_c_bench * aggregate_c_bench);

% compute aggregate welfare in the economy:
aggregate_v_bench = 0;
for i = 1:NL
    for j = 1:NA
        aggregate_v_bench = aggregate_v_bench + V_benchmark(i,j)*psi(i,j);
    end
end

actual_r_bench = alpha * (K_bench / N_bench)^(alpha - 1) - delta;

% results:
K_bench
N_bench
G_bench
Y_bench
aggregate_c_bench
aggregate_v_bench
tau_c_bench

%%%%%%% results:
% ans =
% 
%     0.0529    0.0529
% 
% 
% dist_r =
% 
%    1.7347e-05
% 
% 
% K_bench =
% 
%     7.4287
% 
% 
% N_bench =
% 
%     1.1838
% 
% 
% G_bench =
% 
%     1.1192
% 
% 
% Y_bench =
% 
%     2.4679
% 
% 
% aggregate_c_bench =
% 
%     7.4287
% 
% 
% aggregate_v_bench =
% 
%   -33.0068
% 
% 
% tau_c_bench =
% 
%     0.0750



