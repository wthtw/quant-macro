% the following is our code for solving for an equilibrium for the modified
% Aiyagari model with infinitely-lived households, labour endowment shocks
% l(t) following an AR(1) process, where HH's solve both the
% consumption-savings problem and labour supply problem. HH's can save via
% capital K at interest rate r (endogenously determined) but cannot borrow.
% HH's have convex preferences and the government imposes labour and
% capital income tax as well as consumption tax. The government uses tax
% revenues to finance an exogenous level of G (government spending) and a
% representative firm maximizes their profits given aggregate capital K
% aggregate labour N, the wage rate w and depreciation delta.

clear

% below we list our parameters:
z = 1;  % total factor productivity
alpha = 0.4;    % production function parameter (share of production due to capital)
delta = 0.08;   % proportion of capital saved today for the next period
beta  = 0.96;   % discount factor
rho   = 0.90; % parameters of labor endowment shock process l(t)
              % rho being the autocorrelation coefficient for the AR(1)
              % process:  log(l(t+1)) = rho * log(l(t)) + epsilon(t)
sigma = 0.20; % epsilon(t) is normally distributed with mean zero and
              % standard deviation sigma.
              
tau_y_bench = 0.3;  % labour and capital income tax rate for benchmark
tau_c_bench = 0.075;    % consumption tax rate for benchmark
lambda = 2; % utility function parameter for HH preferences
mu = 0.10;  % parameter used for determining equilibrium interest rate

NA = 400;   % number of intervals in A grid-space, for assets (analogous to K).
NL  = 5; % number of "l" states, for labour efficiency endowment (analgous to Z).

% here we form the grid space for A:
al = 0;
ah = 25;
A = [al : (ah - al )/(NA - 1) : ah];

% here we find a Markov chain which approximates the AR(1) process
% log(l(t+1)) = rho * log(l(t)) + epsilon(t), with mean 0 and s.d. sigma.
% We store this in the transition matrix "pil". These represent the labour
% endowment shocks which the households face.

m = 3; % the bounds of the grid is chosen to be -m*sigma and m*sigma
ll(1) = -m*sigma/(sqrt(1 - rho^2));
ll(NL)=  m*sigma/(sqrt(1 - rho^2));
d = (ll(NL) - ll(1))/(NL - 1); % distance between grids

% here we form the grid space for the labour endowment shocks "l":
for i = 2:NL - 1
    ll(i) = ll(i - 1) + d;
end
l = exp(ll);

% here we form the transition probabilities matrix:
for i = 1:NL
    pil(i,1)  =     normcdf((ll(1) +d/2-rho*ll(i))/sigma,0,1);    
    pil(i,NL) = 1 - normcdf((ll(NL)-d/2-rho*ll(i))/sigma,0,1);
    for j=2:(NL - 1)
        pil(i,j) = normcdf((ll(j)+d/2-rho*ll(i))/sigma,0,1)  ...
                 - normcdf((ll(j)-d/2-rho*ll(i))/sigma,0,1);
    end
end
pil

% here we form the stationary distribution and the value function:
PIL  = [1 0 0 0 0];
PIL  = PIL* pil^1000; % stationary distribution of agents over "l"
LBAR = PIL * l';      % mean labor endowment

% we make an initial guess for the value function:
V_benchmark(1:NL,1:NA) = 0;

% we make an initial guess for the stationary distribution:
psi(1:NL,1:NA) = 0;
psi(1,1) = 1;

% we make an initial guess for interest rate:
dist_r = 1;
r = 1/beta - 1 - 0.001;
r = 0.0379;

% below we solve for the equilibrium using value function iteration:
while dist_r > 0.0001
    R = (1 + (r * (1 - tau_y_bench)));
    KBAR = LBAR * ((r+delta)/(z*alpha)) ^(1/(alpha - 1));
    w = (1 - alpha) * (KBAR/LBAR)^alpha;
  
    % solve the consumer's problem with value function iteration
    dist_V = 1; % this is our starting threshhold for the error
    while dist_V > 0.0000001    % we specify a convergence criterion
        for i = 1:NL
            for j = 1:NA
                n = (1/(1 + lambda)) - ((1/(1 + lambda)) * ((R * A(j) - A)/((1 - tau_y_bench) * w * l(i))));
                n_star = min(max(n, 0), 1);
                C = (R /(1 + tau_c_bench))*A(j) + (((1 - tau_y_bench)/(1 + tau_c_bench))* w * l(i) * n_star) - ((1/(1 + tau_c_bench)) * A);
                U = log(C .* (C > 0)) + (lambda * log(1 - n_star)) + beta * pil(i,:) * V_benchmark;
                % note: do not replace n_star above with n, you will get
                % NaN errors.
                [TV(i,j), indxg(i,j)] = max(U);
            end
        end
        dist_V = max(max(abs(V_benchmark-TV)));
        V_benchmark = TV;
    end
    % consumer's problem has been solved.
    
    % here we find the stationary distribution of the agent:
    dist_psi = 1;   % this is our starting threshold for the error
    while dist_psi > 0.0000000001   % we specify a convergence criterion
        for ip = 1:NL   % loop over lprime:
            for jp = 1:NA    % loop over aprime:
               temp_sum = 0;
               for i = 1:NL   % we sum over l:
                   for j = 1:NA   % we sum over a:
                       temp_sum = temp_sum + psi(i,j) * pil(i,ip) * (indxg(i,j) == jp);
                   end
               end
               Tpsi(ip,jp) = temp_sum;
            end
        end
        dist_psi = max(max(abs(psi-Tpsi))) ;
        psi = Tpsi;
    end
    % psi (the stationary distribution for the economy) found.
    
    KBARn = 0;
    LBARn = 0;
    for i = 1:NL  % l
        for j = 1:NA  % a
            KBARn = KBARn + A(j) * psi(i,j);
            LBARn = LBARn + l(i) * psi(i,j);
        end
    end
    % computed K and L (aggregate capital and labour respectively).
    
    % adjust the interest rate until we achieve equilibrium:
    rn = z*alpha*(KBARn/LBAR)^(alpha - 1) - delta;
    [r rn];
    dist_r = abs(r - rn);
    
    % update the interest rate value:
    r = (1 - mu) * r + mu * rn;
    % computed r'
end

% here we specify the HH's optimal decision rules that we've found:
optimal_n = n_star(indxg);  % optimal labour supply decision.
optimal_c = A(indxg);   % optimal consumption-savings decision.

% we specify the equilibrium levels of aggregate capital and labor,
% respectively, noting that both markets clear.
K_bench = KBARn;
N_bench = LBARn;

% here we specify the aggregate income in the benchmark economy:
Y_bench = K_bench.^(alpha)*(N_bench.^(1 - alpha));

% here we compute the aggregate consumption spending in the benchmark
% economy:
aggregate_c_bench = 0;
for i = 1:NL
    for j = 1:NA
        aggregate_c_bench = aggregate_c_bench + optimal_c(i,j)*psi(i,j);
    end
end
% notice that aggregate_c == K, meaning that we have equilibrium in the
% goods market.

% here we compute the aggregate government spending in the benchmark
% economy, where the government has a balanced budget:
G_bench = (tau_y_bench * (Y_bench - delta * K_bench)) + (tau_c_bench * aggregate_c_bench);

% here we compute aggregate welfare in the economy:
aggregate_v_bench = 0;
for i = 1:NL
    for j = 1:NA
        aggregate_v_bench = aggregate_v_bench + V_benchmark(i,j)*psi(i,j);
    end
end

% here we compute r for the benchmark economy using the equilibrium 
% identity. note that this should be equal to the r found above.
actual_r_bench = alpha * (K_bench / N_bench)^(alpha - 1) - delta;
computed_r_bench = r;
computed_w_bench = w;

% results:
K_bench
N_bench
G_bench
Y_bench
aggregate_c_bench
aggregate_v_bench
tau_c_bench


%%%%%%% expected results:
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



