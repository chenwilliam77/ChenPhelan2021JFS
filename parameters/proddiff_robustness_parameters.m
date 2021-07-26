% This script initializes parameters
%
% Written by William Chen, Mar. 2019

% Add paths
initpaths;

% Economic parameters

% These following commented out parameters use gammaA = gammaB = 1.1, but
% generally should deliver the desired targets no matter the level (but
% does assume symmetric gamma)
% s.rho = .04; s.a = 398; s.a_= 0.85 * s.a; s.kappa = 398; s.lambda = 3; % DO NOT CHANGE
% s.rho = .04; s.a = 1; s.a_= 0.85 * s.a; s.lambda = 42; s.kappa = 5100; % DO NOT CHANGE
% s.rho = .04; s.a = .15; s.a_= 0.85 * s.a; s.lambda = 110; s.kappa = 9100; % DO NOT CHANGE

% For asymmetric gamma
% s.rho = .04; s.a = 1; s.a_= 0.85 * s.a; s.lambda = 50; s.kappa = 5100; % DO NOT CHANGE
s.rho = .04; s.a = 1; s.a_= 0.85 * s.a; s.lambda = 57; s.kappa = 5200; % for asymmetric gamma

s.delta = .05;
s.sigA = .02; s.sigB = .02;
s.gammaA = 1.05; s.gammaB = 2;
s.alpha = 0.5; s.psi = 1;
s.eta_tau = 0;
s.tau = 0;
s.LA = 1e3;
s.LB = 1e3;
% s.LB = 1; % leverage constraint for B

% Settings for numerical solution method
set(groot, 'defaultLineLineWidth', 3);
s.N = 120; s.twist = 1.01; s.start = 1e-3; s.end = 1 - 1e-3;
s.chebnodes = 1; % Use Chebyshev nodes as gridpoints
s.lvg_approx = 0; % Approximate competitive equilibrium from below with sequence of leverage constraints
s.lvg_approx_full = 1; % 0 if permissible to have a binding leverage constraint at the endpoint (generally can find equilibrium without this setting)
s.lvg_approx_iter = 20; % number of iterations for leverage constraint approximation
s.lvg_approx_Linc1 = .5; % initial rate of increase for guessed leverage constraints
s.lvg_approx_Linc2 = .1; % second rate of increase when close to competitive equilibrium
s.lvg_approx_switchLinc = 5e-3; % distance from endpoint at which if leverage constraint binds only there, we lower rate of increase in guessed constraints
s.A_lvg = min(max(2 * 1 / s.LA, 10), 50);
s.B_lvg = 9;
s.LB0 = 26.8; % initial leverage constraint when approximating from below
s.time_dif = .8;
s.static_reltol = 1e-8; s.static_abstol = 1e-8;
s.static_tol = 1e-10;
s.Viter = 1000; s.V_L2tol = 1e-4; s.V_Linftol = 1e-6;
s.qmax_wt = 0;
