% This script produces figures used in the paper
% and determines best responses and optimal policy.
%
% Written by William Chen, Jun. 2019

clear; close all;
addpath ../parameters;

% Choose the correct output file to analyze and the filename of output in savefile
finer_welfare_grid = 1;

if finer_welfare_grid
	output_file = '../data/hpcc/baseline/finer_welfare_grid.mat';
	savefile = '../data/hpcc/baseline/analyze_finer_welfare_grid.mat';
else
	output_file = '../data/hpcc/baseline/coarser_welfare_grid.mat';
	savefile = '../data/hpcc/baseline/analyze_coarser_welfare_grid.mat';
end

%% Settings
recompute_welf = 0;
recompute_stat = 0;
recompute_ce = 0;
parallel = 1;
gammaA = 1.05; gammaB = 2;
pareto_wts = linspace(.05, .95, 30)';
% pareto_wts = linspace(.001, .999, 30)';
diagnostics = 1; % decide whether to show diagnostics

% enter in competitive equilibrium, no leverage expVA and expVB
if recompute_ce
    load('../data/gamA1p05_gamB2.mat');
    baseline_parameters;
    s.eta_tau = 0;
    s.tau = 0;
    s.LA = 1e3;
    s.LB = 1e3;
    vAfnct1 = @(x) interp1(grid.eta, welf.vA, x, 'pchip', 'extrap');
    vBfnct1 = @(x) interp1(grid.eta, welf.vB, x, 'pchip', 'extrap');
    [grid, welf, stat] = get_eqm(vAfnct1, vBfnct1, s, 1, 0, 0);
    ce_expVA = welf.expVA;
    ce_expVB = welf.expVB;
else
    ce_expVA = -5.884779356153428;
    ce_expVB = -3.751850063779598;
end

load(output_file);

if recompute_welf == 1
    disp('Recomputing welfare . . .');
    tim = tic();
    welfs = cell(size(grids));
    stats = cell(size(grids));
    if parallel == 1
        parfor i = 1:numel(grids)
            grid = grids{i};
            if recompute_stat == 1
                [stats{i}.pdf, stats{i}.cdf] = ergocalc(grids{i});
            end
            stat = stats{i};
            vA = grid.xi .* grid.eta .* grid.q;
            vB = grid.zeta .* (1 - grid.eta) .* grid.q;
            VA = vA.^(1 - gammaA) ./ (1 - gammaA);
            VB = vB.^(1 - gammaB) ./ (1 - gammaB);
            expVA = integral(@(x) interp1(grid.eta, stat.pdf .* VA, x, 'pchip', 'extrap'), 0, 1);
            expVB = integral(@(x) interp1(grid.eta, stat.pdf .* VB, x, 'pchip', 'extrap'), 0, 1);
            welfs{i}.vA = vA; welfs{i}.vB = vB;
            welfs{i}.VA = VA; welfs{i}.VB = VB;
            welfs{i}.expVA = expVA;
            welfs{i}.expVB = expVB;
        end
    else
        for i = 1:numel(grids)
            if i / 100 == floor(i / 100)
                disp(i / numel(grids));
            end
            grid = grids{i};
            if recompute_stat == 1
                [stats{i}.pdf, stats{i}.cdf] = ergocalc(grids{i});
            end
            stat = stats{i};
            vA = grid.xi .* grid.eta .* grid.q;
            vB = grid.zeta .* (1 - grid.eta) .* grid.q;
            VA = vA.^(1 - gammaA) ./ (1 - gammaA);
            VB = vB.^(1 - gammaB) ./ (1 - gammaB);
            expVA = integral(@(x) interp1(grid.eta, stat.pdf .* VA, x, 'pchip', 'extrap'), 0, 1);
            expVB = integral(@(x) interp1(grid.eta, stat.pdf .* VB, x, 'pchip', 'extrap'), 0, 1);
            welfs{i}.vA = vA; welfs{i}.vB = vB;
            welfs{i}.VA = VA; welfs{i}.VB = VB;
            welfs{i}.expVA = expVA;
            welfs{i}.expVB = expVB;
        end
    end
    toc(tim);
end

% Check pareto weights to have .5
if sum(pareto_wts == .5) == 0
    pareto_wts = sort([pareto_wts; .5]);
end
output_pareto_wts = [pareto_wts(10); .5; pareto_wts(20)];

%% Get out expected welfare
expVA = zeros(size(welfs));
expVB = expVA;
V_converge = expVA;
for i = 1:numel(welfs)
    expVA(i) = welfs{i}.expVA;
    expVB(i) = welfs{i}.expVB;
    if diagnostics == 1
        V_converge(i) = welfs{i}.V_converge;
    end
end

%% Best responses
disp('Computing best responses . . .');
br_A.welf = zeros(size(welfs, 1), 1);
br_A.pol  = zeros(size(welfs, 1), 1);
br_B.welf = zeros(size(welfs, 2), 1);
br_B.pol  = zeros(size(welfs, 2), 1);
for i = 1:size(welfs, 2)
    ind = find(max(max(expVA(:,i))) == expVA(:,i));
    br_A.welf(i) = expVA(ind,i);
    br_A.pol(i)  = LA_grid(ind,i);
    figure(1);
    plot(LA_grid(:,1), expVA(:,i)); hold on
end
for i = 1:size(welfs, 1)
    ind = find(max(max(expVB(i,:))) == expVB(i,:));
    br_B.welf(i) = expVB(i,ind);
    br_B.pol(i)  = LB_grid(i,ind);
    figure(2);
    plot(LB_grid(1,:), expVB(i,:)); hold on
end

%% Optimal Policies
disp('Computing optimal policies . . .');
expVAf = @(x) obj_VA(x, LA_grid, LB_grid, expVA, 'spline');
expVBf = @(x) obj_VB(x, LA_grid, LB_grid, expVB, 'spline');
opt_flags = zeros(size(pareto_wts));
opt_polA = opt_flags;
opt_polB = opt_flags;
opt_welfA = opt_flags;
opt_welfB = opt_flags;
if parallel == 1
    parfor i = 1:numel(pareto_wts)
        wt = pareto_wts(i);
        [pols, ~, opt_flags(i)] = fminsearch(@(x) -(wt * expVAf(x) + (1 - wt) * expVBf(x)), [.1; .1]);
        opt_polA(i) = pols(1);
        opt_polB(i) = pols(2);
        opt_welfA(i) = expVAf(pols);
        opt_welfB(i) = expVBf(pols);
    end
else
    for i = 1:numel(pareto_wts)
        wt = pareto_wts(i);
        [pols, ~, opt_flags(i)] = fminsearch(@(x) -(wt * expVAf(x) + (1 - wt) * expVBf(x)), [.1; .1]);
        opt_polA(i) = pols(1);
        opt_polB(i) = pols(2);
        opt_welfA(i) = expVAf(pols);
        opt_welfB(i) = expVBf(pols);
    end
end

clear grids welfs stats
save(savefile);
disp('Done');

function out = obj_VA(x, LA_grid, LB_grid, expVA, method)
    if x(1) > max(LA_grid(:,1)) || x(2) > max(LB_grid(1,:))
        out = -1e10;
    elseif x(1) < 0 || x(2) < 0
        out = -1e10;
    else
        out = interp2(LA_grid(:,1), LB_grid(1,:)', expVA', x(1), x(2), method);
    end
end
function out = obj_VB(x, LA_grid, LB_grid, expVB, method)
    if x(1) > max(LA_grid(:,1)) || x(2) > max(LB_grid(1,:))
        out = -1e10;
    elseif x(1) < 0 || x(2) < 0
        out = -1e10;
    else
        out = interp2(LA_grid(:,1), LB_grid(1,:)', expVB', x(1), x(2), method);
    end
end
