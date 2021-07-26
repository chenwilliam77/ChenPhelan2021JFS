% This script produces tables and other calculations used in the paper
% for the CES extension
%
% Written by William Chen, May 2021

clear; close all;
addpath parameters;
initpaths;

save_plots = false;

for sval = [2; 3];
    disp(['Producing results for s = ', num2str(sval)]);

	% File names of the output of welfare computations
	% welf_pol_file1 holds the calculation of welfare on a finer grid
	% of leverage constraints between [0, 1] to make Figures 4 and 5.
	% welf_pol_file2 holds the calculation of welfare on a coarser grid
	% that coveres a larger domain since we want to get a sense
	% of the welfare surface for leverage constraints over
	% (LA, LB) in [0, 20] x [0, 20].
	welf_pol_file1 = ['data/hpcc/ces/ces_s', num2str(sval), '/analyze_welfare_grid.mat'];
	welf_pol_file2 = welf_pol_file1;
	ce_eqm_file    = ['data/hpcc/ces/gamA1p05_gamB2_', num2str(sval), '.mat'];
	set(groot, 'defaultLineLineWidth', 3);
	GAMA = 1.05; GAMB = 2;
	MAXOPTLA = 17; % ymax for plot of optimal LA
	fn = 100; % figure number initial location

	% Welfare computations
	load(ce_eqm_file)
	CE_expVA = welf.expVA;
	CE_expVB = welf.expVB;
	load(welf_pol_file1);
	figure(5);
	plot(LA_grid([1:4,8:size(expVA,1)],1), expVA([1:4,8:size(expVA,1)],1)); hold on
	plot(LA_grid(:,1), expVA(:,find(LA_grid(:,1) >= 0.5 - 1e-3, 1))); hold on
	plot(LA_grid(:,1), expVA(:,find(LA_grid(:,1) == 1.0, 1))); hold off
	xlabel('L_A');
	ylabel('E[V_A]');
	legend({'L_E = 0', 'L_E = 0.5', 'L_E = 1'}, 'FontSize', 16);
	savefig = figure(5);
    if save_plots
	    saveas(savefig, ['figures/ces/ces_s', num2str(sval), '/zfigWelfComps2bGrid'], 'epsc');
    end

	figure(6)
	plot(LB_grid(1,[1:4,8:size(expVA,1)]), expVA(1,[1:4,8:size(expVA,1)])); hold on
	plot(LB_grid(1,:), expVA(find(LA_grid(:,1) >= 0.5 - 1e-3, 1),:)); hold on
	plot(LB_grid(1,:), expVA(find(LA_grid(:,1) == 1.0, 1),:)); hold off
	xlabel('L_E');
	ylabel('E[V_A]');
	legend({'L_A = 0', 'L_A = 0.5', 'L_A = 1'}, 'FontSize', 16, 'location', 'southeast');
	savefig = figure(6);
    if save_plots
	    saveas(savefig, ['figures/ces/ces_s', num2str(sval), '/zfigWelfComps1c'], 'epsc');
    end

	load(welf_pol_file2);
	welfs = LA_grid;
	br_A.welf = zeros(size(welfs, 1), 1);
	br_A.pol  = zeros(size(welfs, 1), 1);
	br_B.welf = zeros(size(welfs, 2), 1);
	br_B.pol  = zeros(size(welfs, 2), 1);

	% Construct the best responses for country A, holding fixed
	% country B's policy, from output of welf_pol_file2
	for i = 1:size(welfs, 2)
	    ind = find(max(max(expVA(:,i))) == expVA(:,i));
        if ind > 0
            br_A.welf(i) = expVA(ind,i);
            br_A.pol(i)  = LA_grid(ind,i);
        end
	    figure(1);
	    plot(LA_grid(:,1), expVA(:,i)); hold on
	end

	% Construct the best responses for country B, holding fixed
	% country A's policy, from output of welf_pol_file2
	for i = 1:size(welfs, 1)
	    ind = find(max(max(expVB(i,:))) == expVB(i,:));
        if ind > 0
            br_B.welf(i) = expVB(i,ind);
            br_B.pol(i)  = LB_grid(i,ind);
        end
	    figure(2);
	    plot(LB_grid(1,:), expVB(i,:)); hold on
	end

	% These should be empty
	indsA = find(br_A.pol > 0);
	indsB = find(br_B.pol > 0);

	% These lines should not run b/c indsA and indsB are empty.
	for i = 1:numel(indsA)
	    figure(201);
	    plot(LA_grid(1:(find(br_A.pol(indsA(i)) == LA_grid(:,1))),1), expVA(1:(find(br_A.pol(indsA(i)) == LA_grid(:,1))),indsA(i))); hold on
	end
	for i = 1:numel(indsB)
	    figure(202);
	    plot(LB_grid(1,1:(find(br_B.pol(indsB(i)) == LB_grid(1,:)))), expVB(indsB(i),1:(find(br_B.pol(indsB(i)) == LB_grid(1,:))))); hold on
	end

	figure(3);
	plot(LB_grid(1,:), br_A.pol);
	xlabel('L_E');
	ylabel('Best Response by Country A');

	figure(4);
	plot(LA_grid(:,1), br_B.pol);
	xlabel('L_A');
	ylabel('Best Response by Country E');

	% Plot the Pareto frontier (given some Pareto weight, we determined
	% the optimal choice of LA and LB when interpolating the
	% welfare surface as a function of LA and LB).
	figure(7);
	plot(opt_welfA, opt_welfB, '-x'); hold on
	plot(opt_welfA(1:floor(numel(pareto_wts) / 2)), opt_welfB(1:floor(numel(pareto_wts) / 2)), '-x'); hold on
	plot(CE_expVA, CE_expVB, 'x', 'color', 'k'); hold on
	plot(expVA(1,1), expVB(1,1), 'o', 'color', 'k'); hold off
	legend({'Pareto Weight >= 1/2', 'Pareto Weight < 1/2', 'Competitive Equilibrium', 'Nash'}, 'FontSize', 12);
	xlim([min([opt_welfA; CE_expVA]) - .1, max(opt_welfA) + .1]);
	ylim([min([opt_welfB; CE_expVB]) - .1, max(opt_welfB) + .1]);
	xlabel('Ex-Ante Welfare of Country A');
	ylabel('Ex-Ante Welfare of Country E');
	title('Pareto Frontier');
	savefig = figure(7);
    if save_plots
	    saveas(savefig, ['figures/ces/ces_s', num2str(sval), '/pareto_frontier'], 'epsc');
    end

	figure(8);
	plot(opt_polA, opt_polB, '-x');
	xlabel('L_A');
	ylabel('L_E');
	title('Optimal Leverage Constraints');

	figure(9);
	plot(pareto_wts(opt_polA < MAXOPTLA), opt_polA(opt_polA < MAXOPTLA));
	xlabel('Pareto Weight');
	ylabel('Optimal L_A');
	savefig = figure(9);
    if save_plots
	    saveas(savefig, ['figures/ces/ces_s', num2str(sval), '/lvgA_pareto_wts'], 'epsc');
    end

	figure(10);
	plot(pareto_wts, opt_polB);
	xlabel('Pareto Weight');
	ylabel('Optimal L_E');
	savefig = figure(10);
    if save_plots
	    saveas(savefig, ['figures/ces/ces_s', num2str(sval), '/lvgB_pareto_wts'], 'epsc');
    end

	for i = 1:numel(output_pareto_wts)
	    valA = num2str(opt_welfA(pareto_wts == output_pareto_wts(i)));
	    valB = num2str(opt_welfB(pareto_wts == output_pareto_wts(i)));
	    polA = num2str(opt_polA(pareto_wts == output_pareto_wts(i)));
	    polB = num2str(opt_polB(pareto_wts == output_pareto_wts(i)));
	    disp(['With a Pareto weight of ' num2str(output_pareto_wts(i)) ' assigned to country A . . .']);
	    disp(['Optimal leverage constraints: (', polA, ', ' polB, ')']);
	    disp(['Ex-ante welfare: (', valA, ', ' valB, ')' newline]);
	end

	s.gammaA = GAMA;
	s.gammaB = GAMB;
	for i = 1:numel(output_pareto_wts)
	    nash2optA = get_consum_equiv(expVA(1,1), opt_welfA(pareto_wts == output_pareto_wts(i)), 'A', s);
	    ce2optA   = get_consum_equiv(CE_expVA, opt_welfA(pareto_wts == output_pareto_wts(i)), 'A', s);
	    nash2optB = get_consum_equiv(expVB(1,1), opt_welfB(pareto_wts == output_pareto_wts(i)), 'B', s);
	    ce2optB   = get_consum_equiv(CE_expVB, opt_welfB(pareto_wts == output_pareto_wts(i)), 'B', s);
	    ce2nashA = get_consum_equiv(CE_expVA, expVA(1,1), 'A', s);
	    ce2nashB   = get_consum_equiv(CE_expVB, expVB(1,1), 'B', s);
	    disp(['With a Pareto weight of ' num2str(output_pareto_wts(i)) ' assigned to country A . . .']);
	    disp(['Consumption equivalent welfare change (pp) from Nash to optimal: (', num2str(100*nash2optA), '%, ' num2str(100*nash2optB), '%)']);
	    disp(['Consumption equivalent welfare change (pp) from CE to optimal: (', num2str(100*ce2optA), '%, ' num2str(100*ce2optB), '%)']);
	    disp(['Consumption equivalent welfare change (pp) from CE to Nash: (', num2str(100*ce2nashA), '%, ' num2str(100*ce2nashB), '%)' newline]);
	end

	% Remove weird points from opt_welfA
    if sval == 10
		opt_welfA = opt_welfA([1:18, 22:30]);
		opt_welfB = opt_welfB([1:18, 22:30]);
		pareto_wts = pareto_wts([1:18, 22:30]);
    end

	% Compute the Pareto weight closest to Nash by pchip interpolation
	pareto_wt_implied_by_A = interp1(opt_welfA, pareto_wts, expVA(1, 1), 'pchip');
	pareto_wt_implied_by_B = interp1(opt_welfB, pareto_wts, expVB(1, 1), 'pchip');
	disp(['The Pareto weight closest to Nash equilibrium is ' num2str( ...
	    (pareto_wt_implied_by_A + pareto_wt_implied_by_B) / 2) '.']);

	% Compute the welfare gain for the EME when the welfare gain for the AE is
	% 0% (relative to CE) by pchip interpolation.
	pareto_wt_A_is_indifferent = interp1(opt_welfA, pareto_wts, CE_expVA, 'pchip');
	implied_expVB = interp1(pareto_wts, opt_welfB, pareto_wt_A_is_indifferent, 'pchip');
	implied_consum_equivB = get_consum_equiv(CE_expVB, implied_expVB, 'B', s);
	disp(['The welfare gain for the EME on the Pareto frontier ', ...
	    'when the AE achieves 0% welfare gain', newline, 'relative to the ', ...
	    'competitive equilibrium is ', num2str(100 * implied_consum_equivB), '%.']);
	disp(['The associated Pareto weight is ', num2str(pareto_wt_A_is_indifferent)]);

	% Compute the welfare gain for the AE when the welfare gain for the EME is
	% 0% (relative to CE) by pchip interpolation.
	pareto_wt_B_is_indifferent = interp1(opt_welfB, pareto_wts, CE_expVB, 'pchip');
	implied_expVA = interp1(pareto_wts, opt_welfA, pareto_wt_B_is_indifferent, 'pchip');
	implied_consum_equivA = get_consum_equiv(CE_expVA, implied_expVA, 'A', s);
	disp(['The welfare gain for the AE on the Pareto frontier ', ...
	    'when the EME achieves 0% welfare gain', newline, 'relative to the ', ...
	    'competitive equilibrium is ', num2str(100 * implied_consum_equivA), '%.']);
	disp(['The associated Pareto weight is ', num2str(pareto_wt_B_is_indifferent)]);

    fprintf '\n\n\n\n\n\n\n';
end
