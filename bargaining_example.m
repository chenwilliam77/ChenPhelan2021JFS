% This script computes the bargaining example in
% Appendix D.1.
%
% Written by William Chen, Jul. 2019

close all; clear;
addpath parameters;
load('data/gamA1p05_gamB2.mat'); % converged solution
vAfnct1 = @(x) interp1(grid.eta, welf.vA, x, 'pchip', 'extrap');
vBfnct1 = @(x) interp1(grid.eta, welf.vB, x, 'pchip', 'extrap');
baseline_parameters;
s.N = 120; s.eta_tau = 0; s.tau = 0; s.gammaA = 1.05; s.gammaB = 2;
s.LA = 0; s.LB = 0;
[~, welf, ~] = get_eqm(vAfnct1, vBfnct1, s, 1, 1, 1);
nol_expVA = welf.expVA;
nol_expVB = welf.expVB;
s.LA = 1e3; s.LB = 1e3;
[~, welf, ~] = get_eqm(vAfnct1, vBfnct1, s, 1, 1, 1);
ce_expVA = welf.expVA; ce_expVB = welf.expVB;
s.eta_tau = .5; s.tau = 4e-3; s.LB = 0;
[grid, welf, ~] = get_eqm(vAfnct1, vBfnct1, s, 1, 1, 1);
interp_pt = .4875;
VA_bargain = interp1(grid.eta, welf.VA, interp_pt, 'linear');
VB_bargain = interp1(grid.eta, welf.VB, interp_pt, 'linear');
disp(['The eta we use is ' num2str(interp_pt)]);
phiA = get_consum_equiv(nol_expVA, VA_bargain, 'A', s);
phiB = get_consum_equiv(nol_expVB, VB_bargain, 'B', s);
disp(['The consumption equivalent change (bp) relative to Nash is (' num2str(phiA*1e4) ...
       ', ' num2str(phiB*1e4) ')']);
phiA_ce = get_consum_equiv(ce_expVA, VA_bargain, 'A', s);
phiB_ce = get_consum_equiv(ce_expVB, VB_bargain, 'B', s);
disp(['The consumption equivalent change (bp) relative to CE is (' num2str(phiA_ce*1e4) ...
       ', ' num2str(phiB_ce*1e4) ')']);
