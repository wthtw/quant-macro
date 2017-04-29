% here we run both benchmark and reform models, finding their equilibrium
% values and comparing them. We also graph the change in aggregate welfare
% changes as a function of asset endowment for each labour endowment state 
% to show the distribution of welfare gains between benchmark and reform
% models. we also graph the value functions for both benchmark and reform
% economies for the sake of completion.

project_3_benchmark_model;
project_3_reform_model;

display('Aggregates for the benchmark economy:')
display('    K:        N:        G:        Y:        agg C:   agg V:     tau_c:')
display('-----------------------------------------------------------------------')
disp([K_bench N_bench G_bench Y_bench aggregate_c_bench aggregate_v_bench tau_c_bench])
display('Equilibrium interest rate and wage rate for the benchmark economy:')
display('    r:        w:')
display('--------------------')
disp([computed_r_bench computed_w_bench])

display('Aggregates for the reform economy:')
display('    K:        N:        G:        Y:        agg C:   agg V:     tau_c:')
display('-----------------------------------------------------------------------')
disp([K_reform N_reform G_reform Y_reform aggregate_c_reform aggregate_v_reform tau_c_reform])
display('Equilibrium interest rate and wage rate for the reform economy:')
display('    r:        w:')
display('--------------------')
disp([computed_r_reform computed_w_reform])

pchangeK = ((K_reform - K_bench) ./ K_bench) .* 100;
pchangeY = ((Y_reform - Y_bench) ./ Y_bench) .* 100;
pchangeaggC = ((aggregate_c_reform - aggregate_c_bench) ./ aggregate_c_bench) .* 100;
pchangeaggV = abs((aggregate_v_reform - aggregate_v_bench) ./ aggregate_v_bench) .* 100;
pchanger = ((computed_r_reform - computed_r_bench) ./ computed_r_bench) .* 100;
pchangew = ((computed_w_reform - computed_w_bench) ./ computed_w_bench) .* 100;
display('% changes between benchmark and reform economies:')
display('   K:         Y:       aggr C:   aggr V:   r:         w:')
display('------------------------------------------------------------')
disp([pchangeK pchangeY pchangeaggC pchangeaggV pchanger pchangew])

for i = 1:NL
    deltaV(i,:) = V_reform(i,:) - V_benchmark(i,:);
    [M, index(i)] = min(deltaV(i,:) .* (deltaV(i,:)));
    indexV(i,:) = A(index(i));
end

clf
axes1 = axes(...
    'FontName', 'Helvetica',...
    'FontSize', 16);
hold on

figure(1)
for i = 1:NL
    plot(A, deltaV(i,:),':','LineWidth', 2)
end
legend('1', '2', '3', '4', '5')
for i = 1:NL
    plot(A(index(i)),deltaV(i,index(i)),'r*')
end
plot([0 25], [0 0], 'k-');  % x-axis
xlabel('A')
ylabel('delta V')
title('Welfare gains for each state (1 to 5)')

figure(2)
hold on
plot(A, optimal_c,':','LineWidth', 2)
plot(A, A,'LineWidth', 1)
xlabel('A')
ylabel('A''')
title('Value functions (benchmark)')

figure(3)
hold on
plot(A, optimal_c_reform,':','LineWidth', 2)
plot(A, A,'LineWidth', 1)
xlabel('A')
ylabel('A''')
title('Value functions (reform)')