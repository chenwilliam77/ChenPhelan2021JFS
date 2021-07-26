% This script computes equilibria used to compare w/benchmark
% cases of linear investment and 1 good, as well as a few other
% use cases, such as comparing competitive eqm to Nash
% when given specific initial condtions.
%
% Written by William Chen, Oct. 2019

load('data/gamA1p05_gamB2.mat'); % converged solution
addpath parameters;
baseline_parameters;
vAfnct1 = @(x) interp1(grid.eta, welf.vA, x, 'pchip', 'extrap');
vBfnct1 = @(x) interp1(grid.eta, welf.vB, x, 'pchip', 'extrap');
LAs = [1e3, 1, .5, .1];
names = cell(4,1);
names{1} = 'ce.mat';
names{2} = 'LA1LBno.mat';
names{3} = 'LAp5LBno.mat';
names{4} = 'LAp1LBno.mat';
for i = 1:4
    s.LA = LAs(i); s.LB = 1e3;
    [grid, welf, stat] = get_eqm(vAfnct1, vBfnct1, s, 1, 0, 0);
    save(['data/', names{i}], 'grid', 'welf', 'stat');
end

% Now compute Nash and a couple other equilibria
s.LA = 0; s.LB = 0;
[grid, welf, stat] = get_eqm(vAfnct1, vBfnct1, s, 1, 0, 0);
save('data/nash.mat', 'grid', 'welf', 'stat');

s.LA = .25; s.LB = .25;
[grid, welf, stat] = get_eqm(vAfnct1, vBfnct1, s, 1, 0, 0);
save('data/LAp25LBp25.mat', 'grid', 'welf', 'stat');

s.LA = .5; s.LB = .5;
[grid, welf, stat] = get_eqm(vAfnct1, vBfnct1, s, 1, 0, 0);
save('data/LAp5LBp5.mat', 'grid', 'welf', 'stat');

s.LA = 1; s.LB = 1;
[grid, welf, stat] = get_eqm(vAfnct1, vBfnct1, s, 1, 0, 0);
save('data/LA1LB1.mat', 'grid', 'welf', 'stat');

s.LA = 2; s.LB = 2;
[grid, welf, stat] = get_eqm(vAfnct1, vBfnct1, s, 1, 0, 0);
save('data/LA2LB2.mat', 'grid', 'welf', 'stat');
