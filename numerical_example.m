% This script provides an example of solving equilibrium
% using this repository's code.
%
% Written by William Chen, Mar. 2019

clear; close all;
close all;
addpath parameters;

% Use saved guess and interpolation, or directly
% guess the value functions?
use_saved_guess = 1;

if saved_guess
	% Load saved guesses and use interpolation
	load('data/gamA1p05_gamB2.mat'); % converged solution
	baseline_parameters;
	vAfnct1 = @(x) interp1(grid.eta, welf.vA, x, 'pchip', 'extrap');
	vBfnct1 = @(x) interp1(grid.eta, welf.vB, x, 'pchip', 'extrap');
    [grid, welf, stat] = get_eqm(vAfnct1, vBfnct1, s, 1, 1, 1);
else
	baseline_parameters;
    [grid, welf, stat] = get_eqm(@(eta) vAfnct(eta, 1, s), @(eta) vBfnct(eta, 1, s), s, 1, 1, 1);
end

display_results(grid, welf, stat, 1);


function out = vAfnct(eta, qweight, s)
    qmin = (s.kappa * sqrt(s.a * s.a_) / 2 + 1) / (s.rho * s.kappa + s.lambda);
    qmax = (s.kappa * s.a / 2 + 1) / (s.rho * s.kappa + s.lambda);
    q = qweight * qmin + (1-qweight) * qmax;
    if s.psi == 1
        out = s.a * q .* eta * s.gammaA / (s.gammaA + s.gammaB) * 2;
    else
        out = s.a * s.rho^(1/(1-s.psi)) * q .* eta * s.gammaA / (s.gammaA + s.gammaB) * 2;
    end
end
function out = vBfnct(eta, qweight, s)
    qmin = (s.kappa * sqrt(s.a * s.a_) / 2 + 1) / (s.rho * s.kappa + s.lambda);
    qmax = (s.kappa * s.a / 2 + 1) / (s.rho * s.kappa + s.lambda);
    q = qweight * qmin + (1-qweight) * qmax;
    if s.psi == 1
        out = s.a * q .* (1-eta) * s.gammaB / (s.gammaA + s.gammaB) * 2;
    else
        out = s.a * s.rho^(1/(1-s.psi)) * q .* (1-eta) * s.gammaB / (s.gammaA + s.gammaB) * 2;
    end
end
