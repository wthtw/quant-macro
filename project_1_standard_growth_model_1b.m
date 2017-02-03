clear

% below is the standard code we used in part 1.a.
% consult that file for documentation. in this file, we're concerned with
% simulating the modeled economy with different levels of initial capital
% stock, and graphing the results.

beta = 0.98;
r = 1/beta - 1;
alpha = 0.36;
z = 1;
delta = 1;

N = 101;
k_ss = (((1/beta) - 1 + delta)* (1/(alpha * z))).^(1 / (alpha - 1));
K = [(0.8 * k_ss):((1.2 * k_ss)-(0.8 * k_ss))/(N - 1):(1.2 * k_ss)]';

C = (1 - delta) * K + (K.^alpha);
V(:) = log(C .*(C > 0));

err_crit = 1;
while (err_crit > 0.00001)
    for i = 1:N
        C = (1 - delta) * K(i) + ((K(i)).^alpha) - K;
        U = log(C .*(C > 0)) + beta * V(:);
        [Ustar, jstar] = max(U);
        TV(i) = Ustar;
        indxg(i) = jstar;
    end
    
    err_crit = max(abs(V - TV));
    V = TV;
end

% simulating the economy with  an initial capital stock 
% of k_1 = coeff *k_ss, where coeff is specified below.

coeff = 0.88;     % here we have a choice of 0.88 or 1.2
ksim(1) = (coeff * k_ss);   % our initial capital stock.
[null, loc_k(1)] = min(abs(K - ksim(1)));    % retrieve the location of our
                                            % initial capital stock, by
                                            % finding the minimum of all
                                            % distances from the entries of
                                            % K to the initial capital
                                            % stock.

for t = 1:39        % we specify 40 periods, since the functions converge
                    % rather quickly.
    loc_k(t + 1) = indxg(loc_k(t));     % we  recursively build a 40-entry 
                                        % vector of locations of subsequent
                                        % values of K (starting from our
                                        % initial capital stock, onwards)
                                        % using the values in indxg as
                                        % reference.
end

ksim = K(loc_k)     % similarly to part 1a, we use the index vector we 
                    % built here, loc_k, to specify the numerical solution.
                    % this will be shown in blue.
                    
k_analytic = alpha * beta * (K(loc_k) .^alpha)  % similarly as before, we
                                                % show the analytic
                                                % solution in red.

clf
axes1 = axes(...
    'FontName', 'Helvetica',...
    'FontSize', 18);
hold on
plot(1:40, ksim, 'LineWidth',4)
plot(1:40, k_analytic, ':', 'LineWidth', 3)
xlabel('t')
ylabel('k(t)')

legend('Numerical', 'Analytical')
%print -dpdf figure1b2.pdf
