% here we do the same calculations as in the benchmark model, but for the
% reform economy. we shall only specify what differs from the benchmark
% in the comments below.

% specify parameters:
z = 1;
alpha = 0.4;
delta = 0.08;
beta  = 0.96;
rho   = 0.90;
sigma = 0.20; 

tau_y_reform = 0;   % here we set the labour and capital income tax rate
                    % for the reform economy as 0.
tau_c_reform = 0.1507;  % here we set the consumption tax for the reform
                        % economy according to the definition:
                        % G_bench / aggregate_c_bench
                        % i.e. 1.1192 / 7.4287
lambda = 2;
mu = 0.10;
G_bench = 1.1192

NA = 400;
NL  = 5;

% form the grid space:
al = 0;
ah = 25;
A = [al : (ah - al )/(NA-1) : ah];


% find a markov chain:
m   = 3; % the bounds of the grid is chosen to be -m*sigma and m*sigma
ll(1) = -m*sigma/(sqrt(1-rho^2));
ll(NL)=  m*sigma/(sqrt(1-rho^2));
d = (ll(NL)-ll(1))/(NL-1); % distance between grids

% form the grid space for the labour endowment shocks "l":
for i = 2:NL - 1
    ll(i) = ll(i - 1) + d;
end
l = exp(ll);

% form the transition probabilities matrix:
for i = 1:NL
    pil(i,1)  =     normcdf((ll(1) +d/2-rho*ll(i))/sigma,0,1);    
    pil(i,NL) = 1 - normcdf((ll(NL)-d/2-rho*ll(i))/sigma,0,1);
    for j=2:(NL - 1)
        pil(i,j) = normcdf((ll(j)+d/2-rho*ll(i))/sigma,0,1)  ...
                 - normcdf((ll(j)-d/2-rho*ll(i))/sigma,0,1);
    end
end
pil

% here we form the stationary distribution and the value function
PIL  = [1 0 0 0 0];
PIL  = PIL* pil^1000;
LBAR = PIL * l';

% make the initial guess for the value function
V_reform(1:NL,1:NA) = 0;

% make an initial guess for the stationary distribution
psi(1:NL,1:NA) = 0;
psi(1,1) = 1;

% make an initial guess for interest rate
dist_r = 1;
r = 1/beta - 1 - 0.001;
r = 0.0379;

% below we solve for the equilibrium using value function iteration:
dist_phi = 1;
while dist_phi > 0.00001
    while dist_r > 0.0001
        R = (1 + (r * (1 - tau_y_reform)));
        KBAR = LBAR * ((r+delta)/(z*alpha)) ^(1/(alpha - 1));
        w = (1 - alpha)*(KBAR/LBAR)^alpha;

        % Solve the consumer's problem with value function iteration
        dist_V = 1; % this is our starting threshhold for the error
        while dist_V > 0.0000001    % we specify a convergence criterion
            for i = 1:NL
                for j = 1:NA 
                    n = (1/(1 + lambda)) - ((1/(1 + lambda))*((R * A(j) - A)/((1 - tau_y_reform) * w * l(i))));
                    n_star = min(max(n, 0), 1);
                    C = (R /(1 + tau_c_reform))*A(j) + (((1 - tau_y_reform)/(1 + tau_c_reform))* w * l(i) * n_star) - ((1/(1 + tau_c_reform)) * A);    % EDITED
                    U = log(C .* (C > 0)) + (lambda * log(1 - n_star)) + beta * pil(i,:) * V_reform;   % EDITED
                    % note: do not replace n_star above with n, you will get
                    % NaN errors.
                    [TV(i,j), indxg(i,j)] = max(U);
                end
            end
            dist_V = max(max(abs(V_reform - TV)));
            V_reform = TV;
        end
        % consumer's problem has been solved.

        % here we find the stationary distribution of agent 
        dist_psi = 1;   % this is our starting threshold for the error
        while dist_psi > 0.0000000001   % we specify a convergence criterion
            for ip = 1:NL   % lprime
                for jp = 1:NA    % aprime 
                   temp_sum=0;
                   for i = 1:NL   % summing over l
                       for j = 1:NA   % summing over a
                           temp_sum = temp_sum + psi(i,j) * pil(i,ip) * (indxg(i,j) == jp);
                       end
                   end
                   Tpsi(ip,jp) = temp_sum;
                end
            end
            dist_psi = max(max(abs(psi - Tpsi))) ;
            psi = Tpsi;
        end
        % psi found.

        KBARn = 0;
        LBARn = 0;
        for i = 1:NL  % l
            for j = 1:NA  % a
                KBARn = KBARn + A(j) * psi(i,j);
                LBARn = LBARn + l(i) * psi(i,j);
            end
        end
        % computed K' and L'

        rn = z * alpha * (KBARn/LBAR)^(alpha-1) - delta;
        [r rn];
        dist_r = abs(r - rn);

        % update interest rate
        r = (1 - mu) * r + mu * rn;

        % computed r'
    end
    
    optimal_c_reform = A(indxg);
    aggregate_c_reform = 0;
    for i = 1:NL
        for j = 1:NA
            aggregate_c_reform = aggregate_c_reform + optimal_c_reform(i,j)*psi(i,j);
        end
    end
    
    K_reform = 0;
    for i = 1:NL
        for j = 1:NA
            K_reform = K_reform + A(j)*psi(i,j);
        end
    end

    N_reform = 0;
    for i = 1:NL
        for j = 1:NA
            N_reform = N_reform + l(i) * psi(i,j);
        end
    end
    
    Y_reform = K_reform.^(alpha)*(N_reform.^(1 - alpha));
    G_reform = (tau_y_reform * (Y_reform - delta * K_reform)) + (tau_c_reform * aggregate_c_reform);
    if tau_c_reform * aggregate_c_reform > G_bench
        tau_c_new = G_bench / aggregate_c_reform;
        tau_c_reform = sigma * tau_c_new + (1 - sigma)* tau_c_reform;
    end
    dist_phi = abs(G_bench - G_reform);
end 


optimal_n_reform = n_star(indxg);


% notice that aggregate_c == K, meaning that we have equilibrium in the
% goods market.

% compute aggregate welfare in the economy:
aggregate_v_reform = 0;
for i = 1:NL
    for j = 1:NA
        aggregate_v_reform = aggregate_v_reform + V_reform(i,j) * psi(i,j);
    end
end

% here we compute r for the refirn economy using the equilibrium 
% identity. note that this should be equal to the r found above.
actual_r_reform = alpha * (K_reform / N_reform)^(alpha - 1) - delta;
computed_r_reform = r;
computed_w_reform = w;

% results:
K_reform
N_reform
G_reform
Y_reform
aggregate_c_reform
aggregate_v_reform
tau_c_reform


%%%%%%% expected results:
% K_reform =
% 
%     9.1932
% 
% 
% N_reform =
% 
%     1.1838
% 
% 
% G_reform =
% 
%     1.1192
% 
% 
% Y_reform =
% 
%     2.6875
% 
% 
% aggregate_c_reform =
% 
%     9.1932
% 
% 
% aggregate_v_reform =
% 
%   -26.5039
% 
% 
% tau_c_reform =
% 
%     0.1217



