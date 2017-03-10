clear
beta = 0.9;
r = 1/beta - 1;
w = 0.1;

N = 1001;
K = [0:1/(N-1):1]';
T = 10;

C = (1+r) * K + w;
V(T,:) = log(C .* (C > 0)) ;   % value function, 
                               % use for backwards iteration

%t = T-1;    % solve for T - 1

for t=T-1:-1:1
    for i=1:N

        C = (1+r) * K(i) + w - K;

        U = log( C .* (C>0) ) + beta * V(t+1,:)' ;  % beta * value f'n 
                                                    % tomorrow
        [Ustar jstar] = max(U);
        [K(i) K(jstar)];
        V(t,i) = Ustar;     % value function
        indxg(t,i) = jstar;     % decision rule

    end
    KPRIME(t,:) =  K(indxg(t,:));
end

clf
axes1 = axes(...
    'FontName', 'Times New Roman',...
    'FontSize', 24);
hold on

plot(K, KPRIME(t,:))
xlabel('k(t)')
ylabel('k(t+1)')

% simulate a household's savings path

ksim(1) = 0.72;    % first period capital stock ("k simulated")

[xyz loc_k(1)] = min(abs(K - ksim(1)))   % location on grid space 
                                         % closest to 0.7
                                         
for t=1:T-1
    loc_k(t + 1) = indxg(t, loc_k(t));
end    

%ksim = K(loc_k)
% wealth goes down monotonically over life cycle
% input: ksim = K(loc_k), ksim
% 
% clf
% axes1 = axes(...
%     'FontName', 'Times New Roman',...
%     'FontSize', 24);
% hold on
% plot(1:T, ksim, 'LineWidth',4)
% xlabel('t')
% ylabel('k(t)')