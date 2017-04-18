
% Let k(t) denote the equilibrium asset level in period t. 
% compute the evolution of capital stock starting from some arbitrary
% capital stock and shock level

T =70000; % max period of simulations
TP=10000; % transition period
shockz   = rand(T,1);
indxk(1) = 1;
indxz(1) = 1;
k(1)     = K(indxk(1));
for t=1:T
    indxk(t+1) = indxg(indxz(t) , indxk(t)); 
    k(t+1)     = K(indxk(t+1));
    invest(t)  = k(t+1) - (1-delta)*k(t); 
    y(t)       = z(indxz(t)) * (k(t)^alpha);
    c(t)       = (1-delta)*k(t) + z(indxz(t)) * (k(t)^alpha) -k(t+1);
    % generate next period's shock
    indxz(t+1) = NZ + 1 - sum(shockz(t) <= cdfz(indxz(t),:));
end
zseq       = z(indxz);

% compute the percentage deviation around the long-run trend for the shock,
% output, consumption, and investment
dt_z(TP:T)=(zseq(TP:T)/mean(zseq(TP:T))-1);
dt_y(TP:T)=(y(TP:T)/mean(y(TP:T))-1);
dt_c(TP:T)=(c(TP:T)/mean(c(TP:T))-1);
dt_invest(TP:T)=(invest(TP:T)/mean(invest(TP:T))-1);

% Compute standard deviations of these percentage deviations
display('------------------------------')
display('   ')
display('    std(y)    std(c)    std(i)')
display('------------------------------')
disp([std(dt_y(TP:T)) std(dt_c(TP:T)) std(dt_invest(TP:T))])
display('------------------------------')
% Compute correlations of these percentage deviations with deviations in
% output
display('  corr(y,c) corr(y,i) ')
display('------------------------------')
disp([corr(dt_y(TP:T)',dt_c(TP:T)') corr(dt_y(TP:T)',dt_invest(TP:T)') ])
display('------------------------------')

display('corr(y(t),y(t+1))')
disp(corr(dt_y(TP:T)',dt_y(TP-1:T-1)'))

clf
TP2=T-100;
plot(dt_y(TP2:T))
hold on
plot(dt_c(TP2:T),'.-r')%;plot(dt_c(TP2:T),'r')
plot(dt_invest(TP2:T),'g')
plot(dt_z(TP2:T),'k--')
legend('output','consumption','investment','TFP')


TP=1;
i=1;
for T=10000:10000:70000
    dt_y(TP:T)=(y(TP:T)/mean(y(TP:T))-1);

    stdy(i) = std(dt_y(TP:T));
    i=i+1;
end

