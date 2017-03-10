clear

% below we specify our designated parameters: discount factor beta,
% production function parameter alpha, depreciation rate of physical
% capital delta, the technology parameter sigma, the discount rate rho,
% the # of states NZ, and the # of intervals in our grid space, NK.

beta = 0.987;   % discount factor
alpha = 0.4;    % elasticity of output wrt capital
sigma = 0.007;  % markov process parameter (variance)
rho = 0.95;     % markov process parameter (autocorrelation)
delta = 1;  % depreciation rate of capital
NZ = 5  % # of states
NK = 100 % # of intervals in K grid-space.

for i = 1:NZ
    for j = 1:NZ
        if i ~= j
            PI(i,j) = 0.05;
        else
            PI(i,j) = 0.80;
        end
    end
end

% below is our dynamic programming algorithm to solve for the firm's
% optimal savings decision on our grid space. we loop over all possible
% values of the grid space, finding for each the optimal decision rule
% jstar and maximized utility Ustar for that particular k in K.

%K = [0:2/(NK - 1):2]';
Z = [0.9 1 1.1 1.2 1.3]';

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

% generates value functions
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
KPRIME_A = (alpha * beta * exp(Z)*(K .^alpha)')'   % analytic solutions

for i = 1:NZ
    Ca(:,i) = ((1 - delta).*K + exp(Z(i)).*(K.^alpha) - KPRIME(:,i));
    Cb(:,i) = ((1 - delta).*K + exp(Z(i)).*(K.^alpha) - KPRIME_A(:,i));
    Ua(:,i) = log(Ca(:,i));
    Ub(:,i) = log(Cb(:,i));
    Ia(:,i) = KPRIME(:,i) - (1 - delta).*K;
    Ib(:,i) = KPRIME_A(:,i) - (1 - delta).*K;
end

clf
axes1 = axes(...
    'FontName', 'Helvetica',...
    'FontSize', 16);
hold on

plot(K, KPRIME_A,':','LineWidth', 3) % plot of analytic capital growth eqn
plot(K, KPRIME,':','LineWidth', 5) % plot of analytic capital growth eqn
plot(K,K,'LineWidth', 1)   % 45 deg line
legend('Numerical', 'Analytical')
xlabel('K')
ylabel('K''')
title('Value functions')

figure(2)
hold on
plot(K, Ca,':','LineWidth', 3)
plot(K, Cb,':','LineWidth', 5)
plot(K,K,'LineWidth', 1)
legend('Numerical', 'Analytical')
xlabel('K')
ylabel('Consumption (C)')
title('Consumption')

figure(3)
hold on
plot(K, Ia,':','LineWidth', 3)
plot(K, Ib,':','LineWidth', 5)
plot(K,K,'LineWidth', 1)
legend('Numerical', 'Analytical')
xlabel('K')
ylabel('Investment (I)')
title('Investment')
%print -dpdf proj2a.pdf
