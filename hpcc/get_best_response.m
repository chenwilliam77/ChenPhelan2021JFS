% This script computes welfare while varying country A's and country B's
% leverage constraint. To use this script, using plot_best_response.m
% first (or some other script) is required to develop a grid of
% value function guesses.
%
% Written by William Chen, Jun. 2019

%% Set up
clear; close all;
addpath ../parameters;
load('data/plot_best_response.mat');
savefile = 'data/welfare_grid.mat';
eta0 = makegrid(s.start, s.end - s.start, s.N);
[guessA, guessB] = make_guess_grid(eta0, vA, vB);
br_lvgA = lvgA;
br_lvgB = lvgB;
s.tau = 0;
baseline_parameters;

% set up initial leverage constraints and Pareto weights
% L0 = [linspace(0, .01, 11)'; linspace(0.05, 1, 20)'; linspace(1.2, 4, (4 - 1) * 5)'; linspace(4.2, 20, (20 - 4.25) * 5)'];
% L0 = [linspace(0, .01, 6)'; linspace(0.05, .4, 8)'; linspace(.5, 1, 11)'];
L0 = [0; 5e-3; .015; .03; .045; .06; .075; .09; .1; linspace(.15, .35, 3)'; .5; linspace(.7, 1, 4)'];
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
end

%% Compute grids, welfares, and statistics
parfor i = 1:numel(LA_grid)
    try
    [grids{i}, welfs{i}, stats{i}] = get_eqm(vf_grid{i}.vAf, vf_grid{i}.vBf, ss{i}, 1, 0, 0);
    catch
        disp(['Error with (LA,LB) pair: (' num2str(LA_grid(i)) ', ' num2str(LB_grid(i)) ')']);
    end
end
clear vf_grid;
save(savefile, 'grids', 'welfs', 'LA_grid', 'LB_grid', 'L0', '-v7.3'); % save data
disp('Done');

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
