% This script initializes parameters
%
% Written by William Chen, Mar. 2019

% Add paths
initpaths;

% Economic parameters
s.rho = .04; s.a = 1; s.a_= 0.85 * s.a; s.lambda = 57; s.kappa = 5200;
s.delta = .05;
s.sigA = .02; s.sigB = .02;
s.gammaA = 1.05; s.gammaB = 2;
s.alpha = 0.5; s.psi = 1;
s.s = 1.1;
s.eta_tau = 0;
s.tau = 0;
s.LA = 1e3;
s.LB = 1e3;
% s.LB = 1; % leverage constraint for B

% Settings for numerical solution method
% Changes for CES relative to baseline: 
% s.N from 120 to 160,
% s.start from 1e-3 to 5e-3 (and s.end also adjusted symmetrically)
% s.V_Linftol from 1e-6 to 1e-5
% Adjustments are needed to obtain convergence 
% when leverage constraints bind
set(groot, 'defaultLineLineWidth', 3);
s.N = 150; s.twist = 1.01; s.start = 1e-2; s.end = 1 - 1e-2;
s.chebnodes = 1; % Use Chebyshev nodes as gridpoints
s.time_dif = .4; % size of implicit time step to update PDE, ranges from zero to one. Higher values => faster but less likely to converge
s.abs_err = true;
s.static_reltol = 1e-8; s.static_abstol = 1e-8; % tolerances for fsolve in static step
s.static_tol = 1e-10;
s.A_lvg = min(max(2 * 1 / s.LA, 10), 50); % guesses for the leverage of each country
s.B_lvg = 9;                              % when their share of aggregate wealth is close to zero
s.Viter = 1000; s.V_L2tol = 1e-4; s.V_Linftol = 1e-6;
s.qmax_wt = 0;
