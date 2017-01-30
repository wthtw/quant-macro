clear
beta = 0.95;     % used to be 0.9, now we are above 45deg line
                % we are accumulating wealth (letting beta = 0.5
r = 1/0.9 - 1;  % previously r = 1/beta - 1, set beta here to remain 0.9
w = 0.1;

N = 1001;
K = [0:1/(N-1):1]';
T = 10;

C = (1 + r) * K + w;
V(:) = log(C .* (C > 0)) ;   % value function, use for backwards iteration


err_crit = 1;
while (err_crit > 0.00001)
    for i=1:N

        C = (1 + r) * K(i) + w - K;

        U = log( C .* (C>0) ) + beta * V(:) ;  % beta * value f'n tomorrow
        [Ustar jstar] = max(U);
        TV(i) = Ustar;     % new value function
                           % solve until we converge
        indxg(i) = jstar;
    end
    
    err_crit = max(abs(V - TV))
    V = TV;
end

KPRIME = K(indxg)
clf
axes1 = axes(...
    'FontName', 'Times New Roman',...
    'FontSize', 24);
hold on

plot(K, KPRIME, 'LineWidth', 4)
plot(K, K)

ksim(1) = 0.72;
[xyz loc_k(1)] = min(abs(K - ksim(1)));

for t = 1:99
    loc_k(t + 1) = indxg(loc_k(t));
end

ksim = K(loc_k);

clf
axes1 = axes(...
    'FontName', 'Times New Roman',...
    'FontSize', 24);
hold on

plot(1:100, ksim, 'LineWidth', 4)
xlabel('t')
ylabel('k(t)')