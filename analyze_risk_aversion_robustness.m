% This script produces results for robustness on risk aversion coefficients.
% Note that this script assumes Nash is (0, 0). You will
% need to alter the scripts in hpcc to check robustness of Nash
% to different risk aversion coefficients.
%
% Written by William Chen, Mar. 2019

% Load saved guesses and use interpolation
clear; close all;
close all;
addpath parameters;
savefile = 'data/risk_aversion_robustness.mat';
do_computation = 0; % set to 0 if you have already computed the results
if do_computation == 1
	load('data/gamA1p05_gamB2.mat'); % converged solution
	risk_aversion_robustness_parameters;
	gammaA_mat = [1.5, 1.4, 1.25, 1.05, 1.05, 1.01, 1.05, 1.05];
	gammaB_mat = [1.5, 1.6, 1.75, 1.9, 2, 2, 2.1, 3];
	vAfnct1 = @(x) interp1(grid.eta, welf.vA, x, 'pchip', 'extrap');
	vBfnct1 = @(x) interp1(grid.eta, welf.vB, x, 'pchip', 'extrap');
	grids_ce = cell(length(gammaA_mat), 1);
	grids_L  = cell(length(gammaA_mat), 1);
	welfs_ce = cell(length(gammaA_mat), 1);
	welfs_L  = cell(length(gammaA_mat), 1);
	for i = 1:length(gammaA_mat)
		s.gammaA = gammaA_mat(i);
		s.gammaB = gammaB_mat(i);

		s.LA = 1e3; s.LB = 1e3;
		[grid, welf, ~] = get_eqm(vAfnct1, vBfnct1, s, 1, 0, 0);
		grids_ce{i} = grid;
		welfs_ce{i} = welf;

		s.LA = 0; s.LB = 0;
		[grid, welf, ~] = get_eqm(vAfnct1, vBfnct1, s, 1, 0, 0);
        grids_L{i} = grid;
		welfs_L{i} = welf;
	end
else
	load(savefile);
end

% Display results
for i = 1:length(gammaA_mat)
	gA = gammaA_mat(i); gB = gammaB_mat(i);
	out_string = ['For risk aversion pair (', num2str(gA), ', ', num2str(gB), '): '];
	s.gammaA = gA;
	s.gammaB = gB;
	phiA = get_consum_equiv(welfs_ce{i}.expVA, welfs_L{i}.expVA, 'A', s);
	phiB = get_consum_equiv(welfs_ce{i}.expVB, welfs_L{i}.expVB, 'B', s);
	out_string = [out_string, num2str(round(phiA * 100, 3)), ', ' ...
				num2str(round(phiB * 100, 3))];
	disp(out_string);

    ss_ind = find(grids_ce{i}.mu_eta < 0, 1);
    eta_interp = (grids_ce{i}.eta(ss_ind) + grids_ce{i}.eta(ss_ind - 1)) / 2;
    debt_output = interp1(grids_ce{i}.eta, grids_ce{i}.debt_gdp, eta_interp, 'linear');
    disp(['Debt/Output: ', num2str(debt_output)]);
end

if do_computation == 1
	save(savefile);
end
