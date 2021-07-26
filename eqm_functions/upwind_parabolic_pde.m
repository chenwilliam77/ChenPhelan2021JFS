function F = upwind_parabolic_pde(X, R, MU, S, G, V, time_dif)
% F = upwind_parabolic_pde(X, R, MU, S, G, V, time_dif)
% where X, R, MU, S, and G are column vectors of the same length
%
% X = [X(1), X(2) ... X(N)]' is the state space (an increasing grid)
% R is the discount rate minus growth
% MU is a drift vector of length N, with MU(1) >= 0 and MU(N) <= 0
% S is a volatility vector of lenght N, with S(1) = S(N) = 0
% G is a payoff flow of length N
%
% we are solving the value function equation 
% R(X)*F(t,X) = G(X) + MU(X)*F_x(t,X) + S(X)*F_xx(t,X)/2 + F_t(t,X) 
% - Main edit is that we replace S^2(X) with S(X) (assume square is done
%   already) so that we can more easily permit multiple Brownian motions
%
% to solve the equation over a small time interval dt with continuation
% value function V, set time_dif = dt/(1 + dt) (or simply time_dif = dt)
%
% to get the stationary solution of this equation, set time_dif = 1: 
% then V has no effect on the solution.
%
% written by Yuliy Sannikov; updated by William Chen, Mar. 2019

if or(MU(1) < 0, MU(end) > 0)
    disp('error: not true that MU(1) >= 0 and MU(N) <= 0');
end
if not(or(S(1) ~= 0, S(end) ~= 0))
    disp('error: not true that S(1) and S(N) both nonzero');
end

N = length(X); 
dX = X(2:N) - X(1:N-1); 

% Perform upwind scheme w/centered difference on diffusion term
S0 = zeros(N,1); S0(2:N-1) = S(2:N-1)./(dX(1:N-2) + dX(2:N-1)); % approx S / (2 * dx): this term is the S/2 coefficient
DU = zeros(N,1); DU(2:N) = - (max(MU(1:N-1),0) + S0(1:N-1))./dX*time_dif; % up diagonal, MU divided by dX, S0 is S / (2 * dx^2)
DD = zeros(N-1,1); DD = - (max(-MU(2:N),0) + S0(2:N))./dX*time_dif; % down diagonal
% observe: S and MU are zero at endpoints, hence S0 zero at endpts too -> 
% boundary conditions for our PDE

D0 = (1 - time_dif)*ones(N,1) + time_dif*R; % diagonal
D0(1:N-1) = D0(1:N-1) - DU(2:N); D0(2:N) = D0(2:N) - DD; % subtract twice b/c centered diff
A = spdiags(D0,0,N,N) + spdiags(DU,1,N,N) + spdiags(DD(1:N-1),-1,N,N);
F = A\(G*time_dif + V*(1 - time_dif)); % solve linear system
end
