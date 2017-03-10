clear

% below we specify our designated parameters: discount factor beta,
% production function parameter alpha, total factor productivity z,
% proportion of capital today saved for the next period, delta, and
% the number of intervals in our grid space, N.

beta = 0.98;
alpha = 0.36;
z = 1;
delta = 1;
N = 101;

% below is the solution for our steady state capital stock in terms of
% the parameters: alpha, beta, delta, and z. It is in a generalized form
% so we can simply change the parameters above to yield a new k_ss.

k_ss = (((1/beta) - 1 + delta)* (1/(alpha * z))).^(1 / (alpha - 1));

% below is our grid-space, of N-many equally spaced points between
% 0.8 * k_ss and 1.2 * k_ss.

K = [(0.8 * k_ss):((1.2 * k_ss)-(0.8 * k_ss))/(N - 1):(1.2 * k_ss)]';

% below we specify our consumption function C and our value function V,
% the latter of which is used for backwards iteration.

C = (1 - delta) * K + (K.^alpha);
V(:) = log(C .*(C > 0));

% below is our dynamic programming algorithm to solve for the firm's
% optimal savings decision on our grid space. we loop over all possible
% values of the grid space, finding for each the optimal decision rule
% jstar and maximized utility Ustar for that particular k in K.

err_crit = 1;       % specify a starting threshhold 
while (err_crit > 0.00001)      % specify a convergence criterion
    for i = 1:N
        C = (1 - delta) * K(i) + ((K(i)).^alpha) - K;
        U = log(C .*(C > 0)) + beta * V(:); % N-vector of value functions
                                            % we find the maximum among
                                            % these functions below.
                                            
        [Ustar, jstar] = max(U);     % assigns the index of the maximizer
                                    % to variable jstar, and maximum value
                                    % to Ustar
                                    
        TV(i) = Ustar;  % the maximized utility Ustar is used for the next
                        % iteration to solve for the previous-period
                        % problem. we store it in the i-th position of 
                        % the N-vector TV.

        
        indxg(i) = jstar;       % set the current decision rule to
                                % the i-th index in the grid space.
                                % we later use the values in indxg
                                % to specify the numerical solution.
    end
    
    err_crit = max(abs(V - TV));    % increment down by the distance
                                    % between consecutive value functions
                                    % (using the supremum norm). we do 
                                    % this until the error satisfies
                                    % the convergence criterion.
                                    
    V = TV;     % after incrementing down, we assign the current N-vector
                % of value-function approximations, TV, as the next
                % N-vector of value-functions V to use in expression U
                % above.
end

KPRIME = K(indxg)   % numerical solution (shown in blue)
                    % i.e. apply the indices of locations of maximizers
                    % found above to the grid K, the largest one
                    % (corresponding to the last one) will be our optimum.

KPRIME_A = alpha * beta * (K .^alpha)   % analytic solution (shown in red)
                                        % this is specified beforehand
                                        % and is used to see how closely
                                        % we approximate the analytic
                                        % solution with the numerical one
                                        % above.
clf
axes1 = axes(...
    'FontName', 'Helvetica',...
    'FontSize', 18);
hold on

plot(K, KPRIME, 'LineWidth', 4) % plot of numerical capital growth eqn
plot(K, KPRIME_A,':','LineWidth', 3) % plot of analytic capital growth eqn
plot(0.1963*ones(size(K)), K,':','LineWidth', 1) % plot of analytic capital growth eqn
plot(K,K)   % 45 deg line
legend('Numerical', 'Analytical')
xlabel('k')
ylabel('k''')
%print -dpdf figure1a.pdf
