% This script produces results for robustness on the productivity parameter a_h.
% Note that this script assumes Nash is (0, 0). You will need to
% alter the scripts in hpcc to check robustness of Nash
% to different productivity parameters.
%
% Written by William Chen, Mar. 2019

% Load saved guesses and use interpolation
clear; close all;
close all;
addpath parameters;
savefile = 'data/proddiff_robustness.mat';
do_computation = 0; % set to 0 if you have already computed the results
if do_computation == 1
	load('data/gamA1p05_gamB2.mat'); % converged solution
	proddiff_robustness_parameters;
	a_mat = [.75, .85, .95];
	vAfnct1 = @(x) interp1(grid.eta, welf.vA, x, 'pchip', 'extrap');
	vBfnct1 = @(x) interp1(grid.eta, welf.vB, x, 'pchip', 'extrap');
	grids_ce = cell(length(a_mat), 1);
	grids_L  = cell(length(a_mat), 1);
	welfs_ce = cell(length(a_mat), 1);
	welfs_L  = cell(length(a_mat), 1);
	for i = 1:length(a_mat)
		s.a_ = a_mat(i);

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
for i = 1:length(a_mat)
	out_string = ['For a_ = ', num2str(a_mat(i)), ': '];
	phiA = get_consum_equiv(welfs_ce{i}.expVA, welfs_L{i}.expVA, 'A', s);
	phiB = get_consum_equiv(welfs_ce{i}.expVB, welfs_L{i}.expVB, 'B', s);
	out_string = [out_string, num2str(round(phiA * 100, 3)), ', ' ...
				num2str(round(phiB * 100, 3))];
	disp(out_string);
end

if do_computation == 1
	save(savefile);
end
