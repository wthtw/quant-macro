clear
beta=0.5;       % play around with this, we began with beta = 0.9
r=1/0.9-1;      % this started as r = 1/beta - 1
w=0.1;

N=1001;
K=[-0.1:1/(N-1):1]';
[N nx] = size(K)

for i=1:N

    C = (1+r)*K(i) + w-K;
    CP = (1 + r) * K + w;
    
    U = log( C .* (C>0) ) + beta * log(CP .* (CP > 0)) ;
    [Ustar jstar] = max(U);
    [K(i) K(jstar)];
    V(i) = Ustar;
    indxg(i) = jstar;

end
KPRIME = K(indxg);

% analytical solution
% as N -> inf, numerical solution -> analytical solution
KPRIME_A = (beta / (1 + beta)) * ( (1 + r) * K +w - w/(beta * (1 + r)));



clf
axes1 = axes(...
    'FontName', 'Times New Roman',...
    'FontSize', 24);
hold on

plot(K, KPRIME, 'LineWidth', 4)
plot(K, KPRIME_A, 'LineWidth', 2)
xlabel('k')
ylabel('kprime')
legend('Numerical', 'Analytical')