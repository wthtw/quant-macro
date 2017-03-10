clear

% below is the standard code we used in part 1.a.
% consult that file for documentation. in this file, we're concerned with
% simulating the modeled economy with different levels of initial capital
% stock, and graphing the results. in this file, we substitute a new value
% of delta, 0.05, from the previous value of 1.

beta = 0.98;
r = 1/beta - 1;
alpha = 0.36;
z = 1;
delta = 0.05;

N = 101;      % default value is N = 101.
              % set N = 11011 for a closer approximation to
              % the steady state value.
              
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

KPRIME = K(indxg);   % numerical solution

% simulating the economy with 
% an initial capital stock of k_1 = 0.88*k_ss:
% below is similar code to part 1.b. for graphing our solutions.
% see the file for 1.b. for step explanations.

coeff = 0.88;     % here we have a choice of 0.88 or 1.2
ksim(1) = (coeff * k_ss);
[null, loc_k(1)] = min(abs(K - ksim(1)));

for t = 1:99
    loc_k(t + 1) = indxg(loc_k(t));
end

ksim = K(loc_k)

clf
axes1 = axes(...
    'FontName', 'Helvetica',...
    'FontSize', 18);
hold on
plot(1:100, ksim, 'LineWidth',4)
xlabel('t')
ylabel('k(t)')
%print -dpdf figure3-2.pdf
