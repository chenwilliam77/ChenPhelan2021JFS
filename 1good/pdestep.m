function [grid, vA, vB, vAt_err, vBt_err] = pdestep(grid, vA0, vB0, s)
% This function performs the pde update step when solving equilibrium.
%
% Written by William Chen, Mar. 2019

% may need to compute sigvAA here

% Compute coefficients
RA = s.gammaA / 2 * ((grid.psiA * s.sigA + grid.sig_vAA).^2 + ...
    (grid.psiB * s.sigB + grid.sig_vAB).^2) + ...
    s.rho * (log(vA0) - log(s.rho * grid.q .* grid.eta)) - ...
    grid.muK - grid.psiA .* grid.sig_vAA * s.sigA - grid.psiB .* grid.sig_vAB * s.sigB;
RB = s.gammaB / 2 * ((grid.psiA * s.sigA + grid.sig_vBA).^2 + ...
    (grid.psiB * s.sigB + grid.sig_vBB).^2) + ...
    s.rho * (log(vB0) - log(s.rho * grid.q .* (1-grid.eta))) - ...
    grid.muK - grid.psiA .* grid.sig_vBA * s.sigA - grid.psiB .* grid.sig_vBB * s.sigB;
% if s.psi == 1
%     RA = RA + s.rho * (log(vA0) - log(s.rho * grid.q .* grid.eta));
%     RB = RB + s.rho * (log(vB0) - log(s.rho * grid.q .* (1-grid.eta)));
% else
%     RA = RA +  1 / (1 - 1/s.psi) .* (s.rho - (vA0 ./ grid.eta ./ grid.q).^(1 - s.psi));
%     RB = RB +  1 / (1 - 1/s.psi) .* (s.rho - (vB0 ./ (1 - grid.eta) ./ grid.q).^(1 - s.psi));
% end
S = grid.eta.^2 .* (grid.sig_etaA.^2 + grid.sig_etaB.^2);
if not(isfield(s, 'time_dif'))
    time_dif = .8;
else
    time_dif = s.time_dif;
end    

% Update VA, VB . . .
% last parameter is dt*rho if dt is small, can be at most 1 (1 corresponds to policy iteration) 
% it is more agressive to set the last parameter closer to 1, but code may not converge
vA = upwind_parabolic_pde(grid.eta, RA, grid.mu_eta .* grid.eta, S, zeros(grid.dim1,1), vA0, time_dif); 
vB = upwind_parabolic_pde(grid.eta, RB, grid.mu_eta .* grid.eta, S, zeros(grid.dim1,1), vB0, time_dif);
end
