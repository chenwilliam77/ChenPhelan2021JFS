% This script analyzes output from hpcc/iter_nash
%
% Written by William Chen, Jul. 2019

clear; close all;
load('data/hpcc/baseline/iter_nash.mat');

for i = 1:numel(LA0)
	disp(['With an initial fixed L_A = ' num2str(LA0(i)) ', leverage constraints iterate to . . .']);
	for j = 1:N_br_iter
		disp(['Iteration ' num2str(j) ': (' num2str(optA{i}(j,1)) ', ' num2str(optA{i}(j,2)), ')']);
	end
	disp(['Final iteration welfare: (' num2str(valA{i}(j,1)) ', ' num2str(valA{i}(j,2)), ')']);
	disp(newline);
end

for i = 1:numel(LB0)
	disp(['With an initial fixed L_B = ' num2str(LB0(i)) ', leverage constraints iterate to . . .']);
	for j = 1:N_br_iter
		disp(['Iteration ' num2str(j) ': (' num2str(optB{i}(j,1)) ', ' num2str(optB{i}(j,2)), ')']);
	end
	disp(['Final iteration welfare: (' num2str(valB{i}(j,1)) ', ' num2str(valB{i}(j,2)), ')']);
	disp(newline);
end
