% This script plots country A's best response, given country B's policy,
% and it plots country B's best response, given country A's policy.
%
% Written by William Chen, Jun. 2019

addpath ../parameters;
load('data/gamA1p05_gamB2.mat');
baseline_parameters;
pareto_weight = 1/2;
lvgA = [5; 4; 3; 2; 1; .9; .6; .5; .3; .2; .1; 0];
lvgB = [5; 4; 3; 2; 1; .9; .6; .5; .3; .2; .1; 0];
[LAs, LBs] = meshgrid(lvgA, lvgB);
LAs = LAs';	 LBs = LBs'; % so varying dim1 is varying LA, varying dim2 is varying LB
VA = cell(size(LAs));
VB = cell(size(LAs));
vA = cell(size(LAs));
vB = cell(size(LAs));
exp_VA = zeros(size(LAs));
exp_VB = zeros(size(LAs));
vA0  = welf.vA;
vB0  = welf.vB;
eta0 = grid.eta;
vAfnct = @(x) interp1(eta0, vA0, x, 'pchip', 'extrap');
vBfnct = @(x) interp1(eta0, vB0, x, 'pchip', 'extrap');

parfor i = 1:numel(LAs)
	  try
	output = get_welfare(s, vAfnct, vBfnct, LAs(i), LBs(i), 1, pareto_weight);
	VA{i}  = output.VA;
	VB{i}  = output.VB;
    vA{i}  = output.vA;
    vB{i}  = output.vB;
	exp_VA(i) = output.expVA;
	exp_VB(i) = output.expVB;
catch
disp(['error with (LA, LB) = (' num2str(LAs(i)) ', ' num2str(LBs(i)) ')']);
end
end

figure(1);
plot(lvgA, exp_VA(:,1)); hold on
plot(lvgA, exp_VA(:,5)); hold on
plot(lvgA, exp_VA(:,10)); hold on
plot(lvgA, exp_VA(:,end - 3)); hold off

figure(2);
plot(lvgB, exp_VB(1,:)); hold on
plot(lvgB, exp_VB(5,:)); hold on
plot(lvgB, exp_VB(10,:)); hold on
plot(lvgB, exp_VB(end - 3,:)); hold off

disp('done');

save('data/plot_best_response.mat');
