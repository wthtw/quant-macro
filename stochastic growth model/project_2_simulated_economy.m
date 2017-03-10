clear

beta = 0.987;   % discount factor
alpha = 0.4;    % elasticity of output wrt capital
sigma = 0.007;  % markov process parameter (variance)
rho = 0.95;     % markov process parameter (autocorrelation)
delta = 0.012;  % depreciation rate of capital

m = 3;
NK = 100 % # of intervals in K grid-space.
NZ = 12  % # of states
sigma_a = sigma/(sqrt(1-(rho^2)));

supremum = m * sigma_a;
infimum = -supremum;
delta2 = (supremum - infimum)/(NZ - 1);
Z = [infimum:delta2:supremum];  % state space of the NZ-state markov chain
zstar = (sum(Z))/NZ;

%%%%% Test Code: setting grid around the steady state capital.
for i = 1:NZ
    k_ss(:,i) = (((1/beta) - 1 + delta)* (1/(alpha * exp(Z(i))))).^(1 / (alpha - 1));
end
kstar = (sum(k_ss))/NZ;
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

% dynamic programming algorithm below.

V = zeros(NK,NZ);
err_crit = 1;       % specify a starting threshhold 
while (err_crit > 0.00000001)      % specify a convergence criterion
    for i = 1:NK
        for j = 1:NZ
            Y(:,j) =  (exp(Z(j)).*((K(i)).^alpha));
            %C(:,j) = ((1 - delta).* K(i)) + (exp(Z(j)).*((K(i)).^alpha)) - K;
            C(:,j) = ((1 - delta).* K(i)) + Y(:,j) - K;
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

KPRIME = K(indxg);   % numerical solutions      
KPRIME_A = (alpha * beta * exp(Z)'*(K .^alpha)')';   % analytic solutions

% Analytic and numerical solutions for consumption, Investment and Output
% for each state 1 to NZ.

for i = 1:NZ
    Ca(:,i) = ((1 - delta).*K + exp(Z(i)).*(K.^alpha) - KPRIME(:,i));
    Cb(:,i) = ((1 - delta).*K + exp(Z(i)).*(K.^alpha) - KPRIME_A(:,i));
    Ua(:,i) = log(Ca(:,i));
    Ub(:,i) = log(Cb(:,i));
    Ia(:,i) = KPRIME(:,i) - (1 - delta).*K;     % numerical
    Ib(:,i) = KPRIME_A(:,i) - (1 - delta).*K;   % analytic
    Ya(:,i) = exp(Z(i)).*(K.^alpha);
end

% output of K(jstar) for NK = 100
Kopt = K(jstar);    % optimal capital for each state (as a vector)
Copt= Ca(jstar)';    % optimal consumption for each state (as a vector)
Iopt = Ia(jstar)';   % optimal investment for each state (as a vector)

T = 1000;
Kt(1) = 0.8*kstar; % Initial K_1, usually kstar
Zt(1) = Z(1);   % Initial Z_1
Zsim(1) = Z(1);
[null, loc_k(1)] = min(abs(K - Kt(1)));
[null, loc_z(1)] = min(abs(Z - Zt(1)));

%%%%% Beginning the simulation

for t = 1:T-1;
    Zsim(t+1) =  rho*Zsim(t) + randn*sigma;
    [null, loc_z(t + 1)] = min(abs(Z - Zsim(t+1))); % find the closest state
                                                    % on our state-grid to
                                                    % the actual AR(1)
                                                    % process above.
                                                    
    index = indxg(:,loc_z(t+1));                    % get the column in indxg
                                                    % corresponding to the
                                                    % state we're in.

    loc_k(t + 1) = index(loc_k(t));                 % get locations below
end

Zt = Z(loc_z);  % retrieves the time-series for technology
Kt = K(loc_k);  % retrieves the time series for capital
Kpt = KPRIME(loc_k);    % retrieves the time series for next-period capital

% code below generates the time series' for consumption, output and
% investment, using the time series' for technology and capital obtained
% above.
for t = 1:T
    Ct(t) = (((1 - delta).*(Kt(t))) + (exp((Zt(t))).*((Kt(t)).^alpha)) - Kpt(t));
    Yt(t) = ((exp(Zt(t)).*(Kt(t)).^alpha));
    It(t) = (Kpt(t) - ((1 - delta) .* Kt(t)));
end

% deviations of variables from their long-run averages.
lrv = 50;   % removes the transition period 
            %(set to remove the first 50 periods)
% below are the time series' with the transition periods removed.
Ylr = Yt(:,(lrv:1:T));
Clr = Ct(:,(lrv:1:T));
Ilr = It(:,(lrv:1:T));
% code below computes the % deviations from long-run trend.
pdevY = ((Ylr - mean(Ylr))/mean(Ylr))*100;
pdevC = ((Clr - mean(Clr))/mean(Clr))*100;
pdevI = ((Ilr - mean(Ilr))/mean(Ilr))*100;
% code belowcomputes the standard deviation of the % deviations 
% from long-run trend
std(pdevY)
std(pdevC)
std(pdevI)

% this code outputs graphs for cross-correlations
% (requires input from command window).
crosscorr(Yt, Yt)
crosscorr(Yt, Ct)
crosscorr(Yt, It)

clf
axes1 = axes(...
    'FontName', 'Helvetica',...
    'FontSize', 16);
hold on

figure(1)
time=(lrv:1:T);
plot(time, pdevY,':','LineWidth', 2)
plot(time, pdevC,':','LineWidth', 1)
plot(time, pdevI,':','LineWidth', 3)
legend('Y(t)', 'C(t)', 'I(t)');
xlabel('Time (t)')
ylabel('% dev from trend')
title('Time-series')
