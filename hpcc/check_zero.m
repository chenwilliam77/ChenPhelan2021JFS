% This script plots country A's best response, given country B's policy,
% and it plots country B's best response, given country A's policy.
%
% Written by William Chen, Jun. 2019

load('../data/gamA1p05_gamB2.mat');
savefile = '../data/hpcc/baseline/check_zero';
addpath parameters;
baseline_parameters;
s.V_Linftol = 1e-7;
pareto_weight = 1/2;
L0 = [0; linspace(1e-3, 1e-2, 10)'; linspace(.02, .1, 9)'];
L0 = [L0, L0];
grids = cell(size(L0));
welfs = cell(size(L0));
stats = cell(size(L0));
ss = cell(size(L0));

for i = 1:size(L0,1)
    s.LA = L0(i,1);
    s.LB = 0;
    ss{i,1} = s;
end

for i = 1:size(L0,1)
    s.LA = 0;
    s.LB = L0(i,1);
    ss{i,2} = s;
end

vA0  = welf.vA;
vB0  = welf.vB;
eta0 = grid.eta;
vAfnct = @(x) interp1(eta0, vA0, x, 'pchip', 'extrap');
vBfnct = @(x) interp1(eta0, vB0, x, 'pchip', 'extrap');


parfor i = 1:size(L0,1)
   try
  [grids{i,1}, welfs{i,1}, stats{i,1}] = get_eqm(vAfnct, vBfnct, ss{i,1}, 1, 0, 0);
   catch
      disp(['Error with (LA, LB) pair = (' num2str(L0(i)) ', 0)']);
   end
end
parfor i = 1:size(L0,1)
   try
  [grids{i,2}, welfs{i,2}, stats{i,2}] = get_eqm(vAfnct, vBfnct, ss{i,2}, 1, 0, 0);
   catch
      disp(['Error with (LA, LB) pair = (0, ' num2str(L0(i)) ')']);
   end
end
   save(savefile, 'grids', 'welfs', 'L0', '-v7.3');
disp('Done');
