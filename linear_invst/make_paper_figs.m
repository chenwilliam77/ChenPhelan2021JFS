% This script illustrates the effects of the terms of trade hedge
% on the long run stationary distribution and equilibrium outcomes.
%
% Written by William Chen, Jun. 2019

clear; close all;
load('data/gamA1p05_gamB2.mat'); % converged solution
savefile = 'data/make_paper_figs';
linear_invst_baseline_parameters;
vAfnct1 = @(x) interp1(grid.eta, welf.vA, x, 'pchip', 'extrap');
vBfnct1 = @(x) interp1(grid.eta, welf.vB, x, 'pchip', 'extrap');

L0 = [1e3; 1; .5; .1];
grids = cell(size(L0)); welfs = grids; stats = grids;
s.LB = 1e3;
for i = 1:length(L0)
    s.LA = L0(i);
    [grids{i}, welfs{i}, stats{i}] = get_eqm(vAfnct1, vBfnct1, s, 1, 1, 1);
end
s.LA = 0;
[grid0, welf0, stat0] = get_eqm(vAfnct1, vBfnct1, s, 1, 1, 1);
loadfiles = cell(4,1);
loadfiles{1} = '../data/ce.mat';
loadfiles{2} = '../data/LA1LBno.mat';
loadfiles{3} = '../data/LAp5LBno.mat';
loadfiles{4} = '../data/LAp1LBno.mat';
gridso = cell(size(loadfiles)); welfso = gridso; statso = gridso; % grids "original"
for i = 1:size(loadfiles,1)
    load(loadfiles{i});
    gridso{i} = grid;
    welfso{i} = welf;
    statso{i} = stat;
end

line_types = {'-', '--', ':', '-x'};
for i = 1:size(grids,1)
    display_results(grids{i}, welfs{i}, stats{i}, line_types{i});
end
figure(4);
legend({'CE', 'L_A = 1', 'L_A = 0.5', 'L_A = 0.1'}, ...
    'Location', 'southwest', 'FontSize', 16);

figure(21);
legend({'CE', 'L_A = 1', 'L_A = 0.5', 'L_A = 0.1'}, ...
    'Location', 'southwest', 'FontSize', 16);

figure(22);
legend({'CE', 'L_A = 1', 'L_A = 0.5', 'L_A = 0.1'}, ...
    'Location', 'southeast', 'FontSize', 16);

figure(23);
legend({'CE', 'L_A = 1', 'L_A = 0.5', 'L_A = 0.1'}, ...
    'Location', 'southwest', 'FontSize', 16);

figure(24);
legend({'CE', 'L_A = 1', 'L_A = 0.5', 'L_A = 0.1'}, ...
    'Location', 'northeast', 'FontSize', 16);

fn = 10000;
for i = 1:2
    figure(fn + 2);
    plot(grids{i}.eta, grids{i}.eta .* sqrt(grids{i}.sig_etaA.^2 + grids{i}.sig_etaB.^2)); hold on
    set(gca, 'ColorOrderIndex', i);
    plot(gridso{i}.eta, gridso{i}.eta .* sqrt(gridso{i}.sig_etaA.^2 + gridso{i}.sig_etaB.^2), '--'); hold on

    figure(fn + 3);
    plot(grids{i}.eta, stats{i}.pdf); hold on
    set(gca, 'ColorOrderIndex', i);
    plot(gridso{i}.eta, statso{i}.pdf, '--'); hold on

    figure(fn + 4);
    plot(grids{i}.eta, grids{i}.eta .* grids{i}.mu_eta); hold on
    set(gca, 'ColorOrderIndex', i);
    plot(gridso{i}.eta, gridso{i}.eta .* gridso{i}.mu_eta, '--'); hold on
end

figure(fn + 2);
xlabel('\eta');
xlim([0, 0.5]);
title('Volatility of \eta');
legend({'CE, Linear', 'CE, Nonlinear', 'L_A = 1, Linear', 'L_A = 1, Nonlinear'}, ...
    'Location', 'northeast');

figure(fn + 3);
xlabel('\eta');
title('Stationary Density');
legend({'CE, Linear', 'CE, Nonlinear', 'L_A = 1, Linear', 'L_A = 1, Nonlinear'});

figure(fn + 4);
xlabel('\eta');
title('Stationary Density');
legend({'CE, Linear', 'CE, Nonlinear', 'L_A = 1, Linear', 'L_A = 1, Nonlinear'});



fn = 20000;
for i = 1:size(gridso,1)
    figure(fn + 1);
    plot(grids{i}.eta, grids{i}.eta .* grids{i}.mu_eta); hold on
    set(gca, 'ColorOrderIndex', i);
    plot(gridso{i}.eta, gridso{i}.eta .* gridso{i}.mu_eta, '--'); hold on

    figure(fn + 2);
    plot(grids{i}.eta, grids{i}.eta .* sqrt(grids{i}.sig_etaA.^2 + grids{i}.sig_etaB.^2)); hold on
    set(gca, 'ColorOrderIndex', i);
    plot(gridso{i}.eta, gridso{i}.eta .* sqrt(gridso{i}.sig_etaA.^2 + gridso{i}.sig_etaB.^2), '--'); hold on

    figure(fn + 3);
    plot(grids{i}.eta, stats{i}.pdf); hold on
    set(gca, 'ColorOrderIndex', i);
    plot(gridso{i}.eta, statso{i}.pdf, '--'); hold on

    figure(fn + 4);
    plot(grids{i}.eta, grids{i}.psiA); hold on
    set(gca, 'ColorOrderIndex', i);
    plot(gridso{i}.eta, gridso{i}.psiA, '--'); hold on
end

figure(fn + 1);
xlabel('\eta');
title('Drift of \eta');
legend({'CE, Linear', 'CE, Nonlinear', 'L_A = 1, Linear', 'L_A = 1, Nonlinear', ...
    'L_A = 0.5, Linear', 'L_A = 0.5, Nonlinear', 'L_A = 0.1, Linear', 'L_A = 0.1, Nonlinear'});

% for i = 1:2
%     figure(fn + 2);
%     set(gca, 'ColorOrderIndex', i);
%     plot(grids{i}.eta, grids{i}.eta .* sqrt(grids{i}.sig_etaA.^2 + grids{i}.sig_etaB.^2)); hold on
%     set(gca, 'ColorOrderIndex', i);
%     plot(gridso{i}.eta, gridso{i}.eta .* sqrt(gridso{i}.sig_etaA.^2 + gridso{i}.sig_etaB.^2), '--'); hold on
% end
figure(fn + 2);
xlabel('\eta');
xlim([0, 0.5]);
title('Volatility of \eta');
% legend({'CE, 1 Good', 'CE, 2 Goods', 'L_A = 1, 1 Good', 'L_A = 1, 2 Goods'});
legend({'CE, Linear', 'CE, Nonlinear', 'L_A = 1, Linear', 'L_A = 1, Nonlinear', ...
    'L_A = 0.5, Linear', 'L_A = 0.5, Nonlinear', 'L_A = 0.1, Linear', 'L_A = 0.1, Nonlinear'}, ...
    'Location', 'northeast');

% for i = 1:2
%     figure(fn + 3);
%     set(gca, 'ColorOrderIndex', i);
%     plot(grids{i}.eta, stats{i}.pdf); hold on
%     set(gca, 'ColorOrderIndex', i);
%     plot(gridso{i}.eta, statso{i}.pdf, '--'); hold on
% end
figure(fn + 3);
xlabel('\eta');
title('Stationary Density');
% legend({'CE, 1 Good', 'CE, 2 Goods', 'L_A = 1, 1 Good', 'L_A = 1, 2 Goods'});
legend({'CE, Linear', 'CE, Nonlinear', 'L_A = 1, Linear', 'L_A = 1, Nonlinear', ...
    'L_A = 0.5, Linear', 'L_A = 0.5, Nonlinear', 'L_A = 0.1, Linear', 'L_A = 0.1, Nonlinear'});

figure(fn + 4);
xlabel('\eta');
title('Share of Capital Held by Country A');
xlim([0,0.5]);
legend({'CE, Linear', 'CE, Nonlinear', 'L_A = 1, Linear', 'L_A = 1, Nonlinear', ...
    'L_A = 0.5, Linear', 'L_A = 0.5, Nonlinear', 'L_A = 0.1, Linear', 'L_A = 0.1, Nonlinear'}, ...
    'Location', 'southeast');

for i = 2:size(gridso,1)
    figure(fn + 5);
    set(gca, 'ColorOrderIndex', i); hold on
    plot(grids{i}.eta, grids{i}.ex_retA); hold on
    set(gca, 'ColorOrderIndex', i);
    plot(gridso{i}.eta, gridso{i}.ex_retA, '--'); hold on

    figure(fn + 6);
    set(gca, 'ColorOrderIndex', i); hold on
    plot(grids{i}.eta, grids{i}.xi .* sqrt(grids{i}.sig_xiA.^2 + grids{i}.sig_xiB.^2)); hold on

    figure(fn + 7);
    set(gca, 'ColorOrderIndex', i); hold on
    plot(gridso{i}.eta, gridso{i}.xi .* sqrt(gridso{i}.sig_xiA.^2 + gridso{i}.sig_xiB.^2)); hold on
end
figure(fn + 5);
xlabel('\eta');
title('Excess Returns for Country A');
legend({'L_A = 1, Linear', 'L_A = 1, Nonlinear', ...
    'L_A = 0.5, Linear', 'L_A = 0.5, Nonlinear', 'L_A = 0.1, Linear', 'L_A = 0.1, Nonlinear'}, ...
    'Location', 'northeast');

figure(fn + 6);
xlabel('\eta');
title('Volatility of \xi');
ylim([0,0.7 * 1e-3]);
legend({'L_A = 1', 'L_A = 0.5', 'L_A = 0.1'}, ...
    'Location', 'northeast', 'FontSize', 16);

figure(fn + 7);
xlabel('\eta');
title('Volatility of \xi');
legend({'L_A = 1', 'L_A = 0.5', 'L_A = 0.1'}, ...
    'Location', 'northeast', 'FontSize', 16);

for i = 1:size(grids,1)
    disp(['Expected welfare of country A: ' num2str(welfs{i}.expVA)]);
    disp(['Expected welfare of country B: ' num2str(welfs{i}.expVB) newline]);
end
disp(['Expected welfare of country A: ' num2str(welf0.expVA)]);
disp(['Expected welfare of country B: ' num2str(welf0.expVB) newline]);

disp('In the 2 good case . . .');
for i = 1:size(gridso,1)
    disp(['Expected welfare of country A: ' num2str(welfso{i}.expVA)]);
    disp(['Expected welfare of country A: ' num2str(welfso{i}.expVB) newline]);
end

save(savefile);

savefig = figure(4);
saveas(savefig, 'figures/tot_LAno_1_p5_p1_LBno', 'epsc');
savefig = figure(21);
saveas(savefig, 'figures/mueta_LAno_1_p5_p1_LBno', 'epsc');
savefig = figure(22);
saveas(savefig, 'figures/voleta_LAno_1_p5_p1_LBno', 'epsc');
savefig = figure(24);
saveas(savefig, 'figures/stat_dist_LAno_1_p5_p1_LBno', 'epsc');
