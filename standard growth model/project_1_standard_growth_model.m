clear
beta = 0.98;
r = 1/beta - 1;
a = 0.36;
z = 1;

N = 101;
k_ss = (1 /(beta * a)).^(1 / (a - 1));
K = [(k_ss * 0.8):((k_ss * 1.2)-(k_ss * 0.8))/(N - 1):(k_ss * 1.2)]';

C = (K.^a);
V(:) = log(C .*(C > 0));

err_crit = 1;
while (err_crit > 0.00001)
    for i = 1:N
        C = ((K(i)).^a) - K;
        U = log(C .*(C > 0)) + beta * V(:);
        [Ustar jstar] = max(U);
        TV(i) = Ustar;
        
        indxg(i) = jstar;
    end
    
    err_crit = max(abs(V - TV));
    V = TV;
end

KPRIME = K(indxg)
clf
axes1 = axes(...
    'FontName', 'Helvetica',...
    'FontSize', 18);
hold on

plot(K, KPRIME, 'LineWidth', 4) % plot of capital growth equation
plot(K,K)   % 45 deg line

