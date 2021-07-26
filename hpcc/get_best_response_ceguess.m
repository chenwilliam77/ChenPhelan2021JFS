% This script plots country A's best response, given country B's policy,
% and it plots country B's best response, given country A's policy.
%
% Written by William Chen, Jun. 2019

addpath ../parameters;
load('data/gamA1p05_gamB2.mat');
savefile = 'data/welfare_grid_ceguess';
baseline_parameters;
pareto_weight = 1/2;
% L0 = [linspace(0, .01, 6)'; linspace(0.05, 0.35, 7)'; linspace(.45, .7, 6)'; linspace(.85, 1, 4)'];
L0 = [0; .0005; .03; .06; .09; .1; .15; linspace(0.25, 0.35, 4)'; .45; .65; .7 ; linspace(.85, 1, 4)'];
[LA_grid, LB_grid] = meshgrid(L0, L0);
LA_grid = LA_grid';
LB_grid = LB_grid'; % so varying dim1 is varying LA, varying dim2 is varying LB
  grids = cell(size(LA_grid));
  welfs = cell(size(LA_grid));
  stats = cell(size(LA_grid));
  ss = cell(size(LA_grid));
for i = 1:numel(LA_grid)
	  s.LA = LA_grid(i);
s.LB = LB_grid(i);
ss{i} = s;
end

vA0  = welf.vA;
vB0  = welf.vB;
eta0 = grid.eta;
vAfnct = @(x) interp1(eta0, vA0, x, 'pchip', 'extrap');
vBfnct = @(x) interp1(eta0, vB0, x, 'pchip', 'extrap');


parfor i = 1:numel(LA_grid)
   try
      [grids{i}, welfs{i}, stats{i}] = get_eqm(vAfnct, vBfnct, ss{i}, 1, 0, 0);
   catch
      disp(['Error with (LA, LB) pair = (' num2str(LA_grid(i)), ', ' num2str(LB_grid(i)) ')']);
   end
end
   save(savefile, 'grids', 'welfs', 'LA_grid', 'LB_grid', 'L0', '-v7.3');
disp('Done');
