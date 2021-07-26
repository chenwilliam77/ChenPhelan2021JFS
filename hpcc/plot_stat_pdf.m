% This script computes welfare while varying country A's and country B's
% leverage constraint. To use this script, using plot_best_response.m
% first (or some other script) is required to develop a grid of
% value function guesses.
%
% Written by William Chen, Jun. 2019

%% Set up
clear; close all;
addpath ../parameters;
load('../data/hpcc/baseline/plot_best_response.mat');
savefile = '../data/hpcc/baseline/plot_stat_pdf.mat';
eta0 = makegrid(s.start, s.end - s.start, s.N);
[guessA, guessB] = make_guess_grid(eta0, vA, vB);
br_lvgA = lvgA;
br_lvgB = lvgB;
baseline_parameters;
s.tau = 0;

% set up initial leverage constraints and Pareto weights
L0 = [0; .1; .5; 1; 30];
[LA_grid, LB_grid] = meshgrid(L0, L0);
LA_grid = LA_grid'; LB_grid = LB_grid';

% make grid of value function guesses
vf_grid = make_vf_grid(guessA, guessB, br_lvgA, br_lvgB, LA_grid, LB_grid);

% make storage grids
grids = cell(size(LA_grid));
welfs = cell(size(LA_grid));
stats = cell(size(LA_grid));
ss    = cell(size(LA_grid));
for i = 1:numel(LA_grid)
    s.LA = LA_grid(i);
    s.LB = LB_grid(i);
    ss{i} = s;
    ss{i}.N = 100;
end

%% Compute grids, welfares, and statistics
parfor i = 1:numel(LA_grid)
    [grids{i}, welfs{i}, stats{i}] = get_eqm(vf_grid{i}.vAf, vf_grid{i}.vBf, ss{i}, 1, 0, 0);
end
clear vf_grid;
save(savefile, 'grids', 'welfs', 'stats', 'ss', 'LA_grid', 'LB_grid', 'L0'); % save data

function out = find_L(LA, LB, lvgA, lvgB)
    out    = zeros(2,1);
	      out(1) = find(lvgA == interp1(lvgA, lvgA, LA, 'nearest', 'extrap'));
	      out(2) = find(lvgB == interp1(lvgB, lvgB, LB, 'nearest', 'extrap'));
end
function vf_grid = make_vf_grid(guessA, guessB, lvgA, lvgB, LA_grid, LB_grid)
    % pre-determine objective functions for fast look up
    getL = @(LA, LB) find_L(LA, LB, lvgA, lvgB);
    vf_grid = cell(size(LA_grid));
    for i = 1:numel(LA_grid)
        inds = getL(LA_grid(i), LB_grid(i));
        vf_grid{i}.vAf = guessA{inds(1), inds(2)};
        vf_grid{i}.vBf = guessB{inds(1), inds(2)};
    end
end
