clear
beta = 0.9;         % discount factor
r = 1/beta - 1;     % interest rate
w = 0.1;

N = 11;     % edit this to make our grid finer
K = [0:1/(N - 1):1]'       % state space for k

% solve for a particular grid point:
% i = 3

k = 0.5;        % this is the value we play around with.
                % if we set k = 0.5, we get negative consumption
                % which isn't possible
                % i.e. C below starts taking negative values
                % we introduce a fix below ***
T = 5       % five periods

% evaluate the obj function for all pts in grid space
% i.e. trying different k' values.

% solve for k = 1 abv but for all pts in K
% generates a column vector of values assumed by U in the grid space
% solve for all values in the K vector (grid space) using a for loop:
for i = 1:N
    C = (1 + r) * K(i) + w - K;     % switch from k to k(i) ***
    U = log(C.* (C > 0)) + beta * log((1 + r) * K + w);     % (C > 0) ***
    
    % graph obj function as a function of k'
    % plot(K, U, '-d')
    
    % find the optimum
    [Ustar jstar] = max(U)
    % K value which maximizes obj function:
    % K(jstar)  % case for k > 0.5, i.e. no negative consumption
    [K(i) K(jstar)]
    V(T - 1, i) = Ustar;        % value function for particular grid pt i
    indxg(T - 1, i) = jstar;    % decision rule, gives location in K vector
                                % of maximizer
                                % we solve the problem for T - 1
end

% locations that maximize obj function:
indxg(T - 1, :);

plot(K, K(indxg(T - 1, :)))     % plots decision rule
                                % will look like a step function.