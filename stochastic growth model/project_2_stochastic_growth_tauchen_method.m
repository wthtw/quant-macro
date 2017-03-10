clear

beta = 0.987;   % discount factor
alpha = 0.4;    % elasticity of output wrt capital
sigma = 0.007;  % markov process parameter (variance)
rho = 0.95;     % markov process parameter (autocorrelation)
delta = 1;  % depreciation rate of capital

m = 3;
NK = 100 % # of intervals in K grid-space.
NZ = 5  % # of states
sigma_a = sigma/(sqrt(1-(rho^2)));

supremum = m * sigma_a;
infimum = -supremum;
delta2 = (supremum - infimum)/(NZ - 1);
Z = [infimum:delta2:supremum];  % state space of the NZ-state markov chain


%%%%% setting grid around the steady state capital.
for i = 1:NZ
    % computing the steady state capital stock for each state.
    k_ss(:,i) = (((1/beta) - 1 + delta)* (1/(alpha * exp(Z(i))))).^(1 / (alpha - 1));
end

kstar = (sum(k_ss))./NZ;    % computes the average steady state capital
kinfimum = 0.8*kstar;
ksupremum = 1.1*kstar;
kstep = (ksupremum - kinfimum)/(NK - 1);
K = [kinfimum:kstep:ksupremum]';
%%%%%


% Here we find a Markov chain which approximates the AR(1) process
% z(t) = rho * z(t-1) + epsilon(t), with mean 0 and s.d. sigma.
% We store this in the transition matrix PI.

for i = 1:NZ
    Fb = ((Z(1) + ((delta2)/2) - rho*Z(i))/sigma);
    PI(i,1) = normcdf(Fb, 0, 1);
end

for i = 1:NZ
    for j = 2:NZ-1
        Fb = ((Z(j) + ((delta2)/2) - rho*Z(i))/sigma);
        Fa = ((Z(j) - ((delta2)/2) - rho*Z(i))/sigma);
        PI(i,j) = normcdf(Fb, 0, 1) - normcdf(Fa, 0, 1);
    end
end

for i = 1:NZ
    Fa = ((Z(NZ) - ((delta2)/2) - rho*Z(i))/sigma);
    PI(i,NZ) = 1 - normcdf(Fa, 0, 1);
end
PI

% below is our dynamic programming algorithm to solve for the firm's
% optimal savings decision on our grid space. we loop over all possible
% values of the grid space, finding for each the optimal decision rule
% jstar and maximized utility Ustar for that particular k in K.

V = zeros(NK,NZ);
err_crit = 1;       % specify a starting threshhold 
while (err_crit > 0.00000001)      % specify a convergence criterion
    for i = 1:NK
        for j = 1:NZ
            C(:,j) = ((1 - delta).* K(i)) + (exp(Z(j)).*((K(i)).^alpha)) - K;
            U(:,j) = log(C(:,j) .* (C(:,j) > 0));
        end
        Util = U + beta * (V*PI);
        [Ustar, jstar] = max(Util);
        TV(i,:) = Ustar;    % column i of TV, Ustar the max row of U(:,:,i)
        indxg(i,:) = jstar; % column i of indxg, jstar a row vector of index
    end
    err_crit = max(max(abs(TV - V)));                    
    V = TV;
end

KPRIME = K(indxg)   % numerical solutions      
KPRIME_A = (alpha * beta * exp(Z)'*(K .^alpha)')'   % analytic solutions

% extra code: simulating evolution equations for capital, utility and
% investment using both numerical and analytic solutions above for
% comparison.

for i = 1:NZ
    Ca(:,i) = ((1 - delta).*K + exp(Z(i)).*(K.^alpha) - KPRIME(:,i));
    Cb(:,i) = ((1 - delta).*K + exp(Z(i)).*(K.^alpha) - KPRIME_A(:,i));
    Ua(:,i) = log(Ca(:,i));
    Ub(:,i) = log(Cb(:,i));
    Ia(:,i) = KPRIME(:,i) - (1 - delta).*K;
    Ib(:,i) = KPRIME_A(:,i) - (1 - delta).*K;
end

% output of K(jstar) for NK = 100
A = K(jstar);   % optimal capital for each state (as a vector)
B = Ca(jstar);  % optimal consumption for each state (as a vector)
D = Ia(jstar);  % optimal investment for each state (as a vector)

clf
axes1 = axes(...
    'FontName', 'Helvetica',...
    'FontSize', 16);
hold on

figure(1)
plot(K, KPRIME_A,':','LineWidth', 2) % plot of analytic capital growth eqn
plot(K, KPRIME,':','LineWidth', 4) % plot of analytic capital growth eqn
for i = 1:NZ
    plot(A(i)*ones(size(K)), K,':', 'LineWidth', 1)
end
plot(K,K,'LineWidth', 1)   % 45 deg line
legend('Numerical', 'Analytical')
xlabel('K')
ylabel('K''')
title('Value functions')

figure(2)
hold on
plot(K, Ca,':','LineWidth', 2)
plot(K, Cb,':','LineWidth', 4)
for i = 1:NZ
    plot(A(i)*ones(size(K)), K,':', 'LineWidth', 1)
end
plot(K,K,'LineWidth', 1)
legend('Numerical', 'Analytical')
xlabel('Capital (k)')
ylabel('Consumption (C)')
title('Consumption')

figure(3)
hold on
plot(K, Ia,':','LineWidth', 2)
plot(K, Ib,':','LineWidth', 4)
for i = 1:NZ
    plot(D(i)*ones(size(K)), K,':', 'LineWidth', 1)
end
plot(K,K,'LineWidth', 1)
legend('Numerical', 'Analytical')
xlabel('K')
ylabel('Investment (I)')
title('Investment')
%print -dpdf proj2b.pdf
