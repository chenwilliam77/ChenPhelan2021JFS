% This script checks that Nash is (0,0) by iteratively finding
% the best response
%
% Written by William Chen, Jul. 2019

clear; close all;
addpath ../parameters;
% load('data/plot_best_response.mat');
load('data/gamA1p05_gamB2.mat');
savefile = 'data/iter_nash_fine_small.mat';
vA0 = welf.vA;
vB0 = welf.vB;
eta0 = grid.eta;
vAf = @(x) interp1(eta0, vA0, x, 'pchip', 'extrap');
vBf = @(x) interp1(eta0, vB0, x, 'pchip', 'extrap');
N_br_iter = 4;
% eta0 = makegrid(s.start, s.end - s.start, s.N);
% [guessA, guessB] = make_guess_grid(eta0, vA, vB);
% br_lvgA = lvgA;
% br_lvgB = lvgB;
parameters;
s.N = 120;

% set up grid of leverage constraints
Lgrid = [0; 1e-3; 2e-3; 3e-3; linspace(0.01, 0.3, 30)'; .5; .7; 1]; % search grid
LA0 = [1; .7; .5; .3; .2; .1; .05; .01]; % initial locations
LB0 = LA0; % start iteration with LB

% gridsA = cell(size(LA0)); % when iteration starts with fixing A's constraint
welfsA = cell(size(LA0));
optA   = cell(size(LA0)); % always L_A first, then L_B
valA   = cell(size(LA0));
% statsA = cell(size(LA0));
ssA    = cell(size(LA0));
% gridsB = cell(size(LB0)); % when iteration starts with fixing B's constraint
welfsB = cell(size(LB0));
optB   = cell(size(LB0)); % always L_A first then L_B
valB   = cell(size(LB0));
% statsB = cell(size(LB0));
ssB    = cell(size(LB0));

% Start with LA0
for i = 1:numel(LA0)
	s.LA = LA0(i);
	welfsA{i} = cell(N_br_iter,1);
	optA{i}   = zeros(N_br_iter,2);
	valA{i}   = zeros(N_br_iter,2);
	for j = 1:N_br_iter
	    welftmp = cell(size(Lgrid));
		if mod(j,2) == 1
			for k = 1:numel(Lgrid)
				% vf_grid = make_vf_grid(guessA, guessB, br_lvgA, br_lvgB, s.LA, Lgrid, 1);
				s.LB = Lgrid(k);
				ssA{k} = s;
			end
		else
			for k = 1:numel(Lgrid)
				s.LA = Lgrid(k);
				ssA{k} = s;
				% vf_grid = make_vf_grid(guessA, guessB, br_lvgA, br_lvgB, s.LB, Lgrid, 0);
			end
		end
		parfor k = 1:numel(Lgrid)
                     try
			% [~, welftmp{k}, ~] = get_eqm(vf_grid{k}.vAf, vf_grid{k}.vBf, ssA{k}, 1, 0, 0);
			[~, welftmp{k}, ~] = get_eqm(vAf, vBf, ssA{k}, 1, 0, 0);
                     catch
                        error(['Error with iteration ' num2str(j) ', Lfree = ' num2str(Lgrid(k)) ]);
                     end
		end
		welfsA{i}{j} = welftmp;
	    expV_tmp = zeros(numel(Lgrid),1);
        if mod(j,2) == 1
            for k = 1:numel(Lgrid)
			    expV_tmp(k) = welfsA{i}{j}{k}.expVB;
            end
        else
            for k = 1:numel(Lgrid)
			    expV_tmp(k) = welfsA{i}{j}{k}.expVA;
            end
        end
		[val,index] = max(reshape(expV_tmp, numel(expV_tmp), 1));
        if mod(j,2) == 1
			s.LB = Lgrid(index);
			valA{i}(j,:) = [welfsA{i}{j}{index}.expVA, val];
        else
			s.LA = Lgrid(index);
			valA{i}(j,:) = [val, welfsA{i}{j}{index}.expVB];
        end
        optA{i}(j,:) = [s.LA, s.LB];
	end
end

% Onto LB0
for i = 1:numel(LB0)
	s.LB = LB0(i);
	welfsB{i} = cell(N_br_iter,1);
	optB{i}   = zeros(N_br_iter,2);
	valB{i}   = zeros(N_br_iter,2);
	for j = 1:N_br_iter
	    welftmp = cell(size(Lgrid));
		if mod(j,2) == 1
			for k = 1:numel(Lgrid)
				s.LA = Lgrid(k);
				ssB{k} = s;
				% vf_grid = make_vf_grid(guessA, guessB, br_lvgA, br_lvgB, s.LB, Lgrid, 0);
			end
		else
			for k = 1:numel(Lgrid)
				s.LB = Lgrid(k);
				ssB{k} = s;
				% vf_grid = make_vf_grid(guessA, guessB, br_lvgA, br_lvgB, s.LA, Lgrid, 1);
			end
		end

		parfor k = 1:numel(Lgrid)
			% [~, welftmp{k}, ~] = get_eqm(vf_grid{k}.vAf, vf_grid{k}.vBf, ssB{k}, 1, 0, 0);
			[~, welftmp{k}, ~] = get_eqm(vAf, vBf, ssB{k}, 1, 0, 0);
		end
		welfsB{i}{j} = welftmp;
	    expV_tmp = zeros(numel(Lgrid),1);
        if mod(j,2) == 1
            for k = 1:numel(Lgrid)
			    expV_tmp(k) = welfsB{i}{j}{k}.expVA;
            end
        else
            for k = 1:numel(Lgrid)
			    expV_tmp(k) = welfsB{i}{j}{k}.expVB;
            end
         end
		[val,index] = max(reshape(expV_tmp, numel(expV_tmp), 1));
        if mod(j,2) == 1
			s.LA = Lgrid(index);
			valB{i}(j,:) = [val, welfsB{i}{j}{index}.expVB];
		else
			s.LB = Lgrid(index);
			valB{i}(j,:) = [welfsB{i}{j}{index}.expVA, val];
        end
		optB{i}(j,:) = [s.LA, s.LB];
	end
end

save(savefile, 'Lgrid', 'LA0', 'LB0', 'optA', 'valA', 'optB', 'valB', 'N_br_iter', '-v7.3');

disp('Done');

function out = find_L(LA, LB, lvgA, lvgB)
    out    = zeros(2,1);
          out(1) = find(lvgA == interp1(lvgA, lvgA, LA, 'nearest', 'extrap'));
          out(2) = find(lvgB == interp1(lvgB, lvgB, LB, 'nearest', 'extrap'));
end
function vf_grid = make_vf_grid(guessA, guessB, lvgA, lvgB, L0, Lgrid, fixA)
    % pre-determine objective functions for fast look up
	if fixA == 1
	    getL = @(LB) find_L(L0, LB, lvgA, lvgB);
	else
	    getL = @(LA) find_L(LA, L0, lvgA, lvgB);
	end
    vf_grid = cell(length(Lgrid))
    for i = 1:numel(Lgrid)
        inds = getL(Lgrid(i));
        vf_grid{i}.vAf = guessA{inds(1), inds(2)};
        vf_grid{i}.vBf = guessB{inds(1), inds(2)};
    end
end
