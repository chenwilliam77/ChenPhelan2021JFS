function display_results(grid, welf, stat, diagnostics, pareto_weight)
% function display_results(grid, welf, stat, diagnostics, pareto_weight)
% This function plots equilibrium objects.
% If diagnostics == 1, then we also plot a variety of
% error diagnostics to verify that convergence is
% achieved.
% In one of the plots, we assign pareto_weight to country A. If
% not provided, we set it to 1/2.
%
% Written by William Chen, Gregory Phelan Jun 2019

if nargin < 4
	diagnostics = 0;
end
if nargin < 5
	pareto_weight = 1/2;
end

figure(1);
plot(grid.eta, grid.q); hold on
xlabel('\eta');
title('Asset Price q');


figure(2);
plot(grid.eta, grid.psiAa); hold on
hFigure2 = findall(0, 'type', 'figure', 'number', 2);
hLines2 = findall(hFigure2, 'type', 'line');
n_lines = size(hLines2, 1);
set(gca, 'ColorOrderIndex', floor(n_lines / 2) + 1);
plot(grid.eta, grid.psiAb, 'LineStyle', '--'); hold on
legend({'\psi_{Aa}', '\psi_{Ab}'}, 'FontSize', 16, 'Location', 'south');
xlabel('\eta')
title('Capital Allocations of Country A');

figure(3);
plot(grid.eta, grid.psiBb); hold on
hFigure3 = findall(0, 'type', 'figure', 'number', 3);
hLines3 = findall(hFigure3, 'type', 'line');
n_lines = size(hLines3, 1);
set(gca, 'ColorOrderIndex', floor(n_lines / 2) + 1);
plot(grid.eta, grid.psiBa, 'LineStyle', '--'); hold on
legend({'\psi_{Bb}', '\psi_{Ba}'}, 'FontSize', 16, 'Location', 'south');
xlabel('\eta');
title('Capital Allocations of Country B');

figure(4);
plot(grid.eta, grid.Pa ./ grid.Pb); hold on
xlabel('\eta');
title('Terms of Trade P_a / P_b');


figure(5);
plot(grid.eta, grid.lvgA); hold on
xlabel('\eta');
title('Leverage of Country A');

figure(6);
plot(grid.eta, grid.lvgB); hold on
xlabel('\eta');
title('Leverage of Country B');

figure(7);
plot(grid.eta, grid.ex_retA); hold on
hFigure6 = findall(0, 'type', 'figure', 'number', 6);
hLines6 = findall(hFigure6, 'type', 'line');
n_lines = size(hLines6, 1);
set(gca, 'ColorOrderIndex', floor(n_lines / 2) + 1);
plot(grid.eta, grid.ex_retB, 'LineStyle', '--'); hold on
xlabel('\eta');
legend({'Country A', 'Country B'}, 'FontSize', 16, 'Location', 'north');
title('Excess Returns');

figure(8);
plot(grid.eta, grid.muK); hold on
xlabel('\eta');
title('Growth Rate of Capital');

figure(9);
plot(grid.eta, grid.rF); hold on
xlabel('\eta');
title('Risk-Free Interest Rate');

figure(21);
plot(grid.eta, grid.eta .* grid.mu_eta); hold on
xlabel('\eta');
title('Drift of \eta');

figure(22);
plot(grid.eta, grid.eta .* sqrt(grid.sig_etaA.^2 + grid.sig_etaB.^2)); hold on
xlabel('\eta');
title('Volatility of \eta');

figure(23);
plot(grid.eta, grid.q .* sqrt(grid.sig_qA.^2 + grid.sig_qB.^2)); hold on
xlabel('\eta');
title('Volatility of q');

figure(24);
plot(grid.eta(2:end-1), stat.pdf(2:end-1)); hold on
xlabel('\eta');
title('Stationary Density');

figure(41);
plot(grid.eta, welf.VA); hold on
xlabel('\eta');
title('Welfare of Country A');

figure(42);
plot(grid.eta, welf.VB); hold on
xlabel('\eta');
title('Welfare of Country B');

figure(43);
plot(grid.eta, pareto_weight * welf.VA + (1 - pareto_weight) * welf.VB); hold on
xlabel('\eta');
title('Average Welfare for Social Planner');

figure(44);
plot(grid.eta, grid.xi); hold on
xlabel('\eta');
title('Marginal Value of Wealth in Country A \xi');

figure(45);
plot(grid.eta, grid.zeta); hold on
xlabel('\eta');
title('Marginal Value of Wealth in Country B \zeta');

figure(46);
plot(grid.eta, grid.xi ./ grid.zeta); hold on
xlabel('\eta');
title('Relative Marginal Value of Wealth \xi/\zeta');


if diagnostics == 1
	fignum = 100;
	figure(fignum + 1);
plot(grid.eta, grid.flag); hold on
xlabel('\eta');
title('fsolve Successes');

	figure(fignum + 2);
plot(grid.eta, welf.vAerr); hold on
xlabel('\eta');
title('Change in v_A from Last Step');

	figure(fignum + 3);
plot(grid.eta, welf.vBerr); hold on
xlabel('\eta');
title('Change in v_B from Last Step');

i_end = size(welf.vAerr_all,1);
for i = 1:size(welf.vAerr_all,2)
    if sum(welf.vAerr_all(:,i)) + sum(welf.vBerr_all(:,i)) == 0
        i_end = i-1;
        break;
    end
end
j = i_end - 3;
for i = j:i_end
    figure(fignum + 4);
    plot(grid.eta, welf.vAerr_all(:,i)); hold on

    figure(fignum + 5);
    plot(grid.eta, welf.vBerr_all(:,i)); hold on
end
figure(fignum + 4); hold off
xlabel('\eta');
title('Change in v_A from Past Few Steps');
figure(fignum + 5); hold off
title('Change in v_B from Past Few Steps');
end
