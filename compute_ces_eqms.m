% This script computes the competitive equilibrium with CES production functions
% using a homotopy approach. The output is then saved to provide
% initial guesses when leverage constraints are imposed in equilibrium.
%
% Written by William Chen, May 2021

% Load saved guesses and use interpolation
clear; close all;
addpath parameters;
load('data/gamA1p05_gamB2.mat'); % converged solution
ces_parameters;
s_vec = [2; 3; 4; 5; 6; 10;];
s_vec_name = {'2', '3', '4', '5', '6', '10'};

for i = 1:numel(s_vec)
    disp(['Solving case of s = ', num2str(s_vec(i))]);
    if i > 1
        clear grid welf stat
        load(['data/hpcc/ces/gamA1p05_gamB2_', s_vec_name{i-1}, '.mat']);
    end
    s.s = s_vec(i);
    vAfnct1 = @(x) interp1(grid.eta, welf.vA, x, 'spline', 'extrap');
    vBfnct1 = @(x) interp1(grid.eta, welf.vB, x, 'spline', 'extrap');
    [grid, welf, stat] = get_eqm(vAfnct1, vBfnct1, s, 1, 1, 1);
    save(['data/hpcc/ces/gamA1p05_gamB2_', s_vec_name{i}, '.mat'], 'grid', 'welf', 'stat');
end
