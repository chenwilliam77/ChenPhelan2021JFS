% This function solves for welfare, given a user input of parameters,
% value function guesses as functions of eta,
% and leverage constraints for countries A and B. The output depends on the user input.
% If type == 1, then we output a struct with value functions, expected welfare, and a Boolean for convergence.
% If type == 2, then we output the sum of value functions.
% If type == 3, then we output the sum of expected welfare.
% If pareto_weight is provided (optional setting), then we sum using these pareto weights.
% Default is 1/2. Weight is placed on country A
%
% Written by William Chen, Jun. 2019
function out = get_welfare(s, vAfnct, vBfnct, LA, LB, type, pareto_weight)
	if nargin == 5
		type = 1;
	end
	if nargin < 7
		pareto_weight = 1/2;
	end

	% solve equilibrium
	s.LA = LA;
	s.LB = LB;

	% Determine out
	if type == 1
		[grid, welf, ~] = get_eqm(vAfnct, vBfnct, s, 1, 0, 0);
		out.eta   = grid.eta;
		out.VA    = welf.VA;
		out.VB    = welf.VB;
		out.expVA = welf.expVA;
		out.expVB = welf.expVB;
        out.V_converge = welf.V_converge;
	elseif type == 2
		[~, welf, ~] = get_eqm(vAfnct, vBfnct, s, 0, 0, 0);
		out = pareto_weight * welf.VA + (1 - pareto_weight) * welf.VB;
	elseif type == 3
		[~, welf, ~] = get_eqm(vAfnct, vBfnct, s, 1, 0, 0);
		out = pareto_weight * welf.expVA + (1 - pareto_weight) * welf.expVB;
	elseif type == 4
		% Country A's expected welfare
		[~, welf, ~] = get_eqm(vAfnct, vBfnct, s, 1, 0, 0);
		out = welf.expVA;
	elseif type == 5
		% Country B's expected welfare
		[~, welf, ~] = get_eqm(vAfnct, vBfnct, s, 1, 0, 0);
		out = welf.expVB;
	else
		error('Incorrect value provided by user for the variable type');
	end
end
