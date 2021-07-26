function [grid,s] = staticstep_fd_ces(grid, vA, vB, s, lastiter)
% [grid, s] = staticstep_fd_ces(grid, vA, vB, s, lasiter)
%
% This function performs the static step to solve for equilibrium,
% given vA and vB, when the production function is CES
%
% Written by William Chen, May 2021

if nargin < 5
    lastiter = 0;
elseif not(lastiter == 1 || lastiter == 0)
    error('Last argument is not 0 or 1');
end

%% Set up

% fsolve settings
options = optimset('display', 'off', 'TolFun', s.static_tol);

% Compute first-order derivatives of vA, vB for equilibrium computations
deta = diff(grid.eta); vAp = diff(vA) ./ deta; vBp = diff(vB) ./ deta;
vAp_f = [vAp; 0]; vAp_b = [0; vAp]; vBp_f = [vBp; 0]; vBp_b = [0; vBp];
vAp_c = [0; (vA(3:end) - vA(1:end-2)) ./ sum(deta(2:end) + deta(1:end-1)); 0];
vBp_c = [0; (vB(3:end) - vB(1:end-2)) ./ sum(deta(2:end) + deta(1:end-1)); 0];

if ~isreal(vA) || ~isreal(vB)
    keyboard;
end

vA_interp    = griddedInterpolant(grid.eta, vA, 'spline', 'spline'); % spline to impose C2 cubic
vB_interp    = griddedInterpolant(grid.eta, vB, 'spline', 'spline');
vAp_f_interp = griddedInterpolant(grid.eta(1:end-1), vAp_f(1:end-1), 'pchip', 'pchip'); % pchip b/c don't need C2
vBp_f_interp = griddedInterpolant(grid.eta(1:end-1), vBp_f(1:end-1), 'pchip', 'pchip'); % and shape-preserving is probably
vAp_b_interp = griddedInterpolant(grid.eta(2:end), vAp_b(2:end), 'pchip', 'pchip');     % best since the derivative
vBp_b_interp = griddedInterpolant(grid.eta(2:end), vBp_b(2:end), 'pchip', 'pchip');     % should have a "regular" shape
vAp_c_interp = griddedInterpolant(grid.eta(2:end-1), vAp_c(2:end-1), 'pchip', 'pchip');
vBp_c_interp = griddedInterpolant(grid.eta(2:end-1), vBp_c(2:end-1), 'pchip', 'pchip');

% Store interpolation functions
v_interp.vA    = vA_interp;
v_interp.vB    = vB_interp;
v_interp.vAp_f = vAp_f_interp;
v_interp.vBp_f = vBp_f_interp;
v_interp.vAp_b = vAp_b_interp;
v_interp.vBp_b = vBp_b_interp;
v_interp.vAp_c = vAp_c_interp;
v_interp.vBp_c = vBp_c_interp;

% Find threshold capital allocation of psiAa or psiBb
% after which specialization stops
s.threshold_cap_Aa = 1 / (1 + (s.a / s.a_ * (1 - s.alpha) / s.alpha)^s.s);
s.threshold_cap_Bb = 1 / (1 + (s.a / s.a_ * s.alpha / (1 - s.alpha))^s.s);

%% Solve from the left first
% Make guess for capital allocation
eta0   = grid.eta(1);
psiAa0 = s.A_lvg * eta0;
psiBa0 = (1 - (1 + (s.a / s.a_ * (1 - s.alpha) / s.alpha)^s.s) * psiAa0) / ...
         (1 + (s.a / s.a_)^(s.s - 1) * ((1 - s.alpha) / s.alpha)^s.s);
psiBb0 = 1 - psiAa0 - psiBa0;
Y = ces_production((s.a * psiAa0 + s.a_ * psiBa0), s.a * psiBb0, s.alpha, s.s);
q0 = (s.kappa * Y + 1) / (s.rho * s.kappa + s.lambda);
if s.psi ~= 1
    q0 = fsolve(@(q) get_q_eis(q, eta0, vA(1), vB(1), Y, s), q0, optimset('display', 'off'));
end

% Derive minimum capital price for left-hand side.
% The choice of psiBa is obtained from formula for psiBa
% when EME stops specializing by setting psiAa = 0
Ymin_psiBa = 1 / (1 + (s.a / s.a_)^(s.s - 1) * ((1 - s.alpha) / s.alpha)^s.s);
Ymin_psiBb = 1 - Ymin_psiBa;
Ymin = ces_production(s.a_ * Ymin_psiBa, s.a * Ymin_psiBb, s.alpha, s.s);
if s.psi ~= 1
    qmin = fsolve(@(q) get_q_eis(q, 0, vA(1), vB(1), Ymin, s), q0, optimset('display', 'off'));
else
    qmin = (s.kappa * Ymin + 1) / (s.rho*s.kappa + s.lambda);
end
s.qminA = qmin; % if psi ~= 1, the qmin may be different when starting from A's side (left)

% Find eta which maximizes q to figure out where to
% split the state space in three regions.
% To the left of eta_maxq, we expect q to increase with eta.
% To the right, we expect q to decrease with eta.
% Region 1 (B does not specialize) should occur left of eta_maxq,
% region 2 (both specialize) should occur to both the left and right,
% and region 3 (A does not specialize) should occur right of eta_maxq.
opt_Ya = s.a * (s.alpha^s.s) / ((1 - s.alpha)^s.s + s.alpha^s.s);
opt_Yb = s.a * (1-s.alpha)^s.s / ((1 - s.alpha)^s.s + s.alpha^s.s);
qmax = (s.kappa * ces_production(opt_Ya, opt_Yb, s.alpha, s.s) + 1) / ...
        (s.rho * s.kappa + s.lambda);
s.qmax = qmax;
if s.gammaA == s.gammaB
    n_eta = ceil(s.N / 2);
else
    initpt = s.alpha * (s.gammaB / (s.gammaA + s.gammaB));
    eta_maxq = fsolve(@(eta) find_argqmax(eta, v_interp, s), initpt, options);
    eta0 = interp1(grid.eta, grid.eta, eta_maxq, 'nearest');
    n_eta = find(grid.eta == eta0); % location of the initial starting eta, which maximizes q
    if eta0 == eta_maxq % if there's a grid point exactly equal to eta_maxq
	    grid.eta(n_eta) = eta0 - (eta0 - grid.eta(n_eta-1)) / 3; % then decrease it by a tiny amount
	    % eta0 = grid.eta(n_eta);
	elseif eta0 > eta_maxq % if eta0 > eta_maxq, decrease n_eta b/c we want
	    n_eta = n_eta - 1; % n_eta to be the index of the last grid point for which
	    % eta0 = grid.eta(n_eta); % q is increasing with eta
    end

    % one last adjustment since the algorithm seems to
    % work better if a greater portion of region 2 is solved
    % starting from the left of the state space,
    % so we take a convex combo of n_eta and s.N / 2,
    % but with more weight on n_eta so it's not too far
    % from the eta for which q is maximized.
    n_eta = ceil((2 / 3) * n_eta + (1 / 3) * s.N / 2);
end

%% Solve using backward differences
% start with region 1, backward diff
fd_vals = [qmin; s.start];
grid.flag = zeros(size(grid.eta));
s.side = 0; % starting before max q
for j = 1:n_eta
    eqc = @(x) eqcond(grid.eta(j), x, fd_vals, v_interp, grid, s, 1);

    % Transform initial guesses
    transformed_q0   = square_root_interval_to_real_line(q0, s.qminA, s.qmax);
    transformed_psi0 = square_root_interval_to_real_line(psiAa0, 0, 1); % for side = 0, psiAa0 = x(2)

    % Solve!
    [out1, ~, flag, ~] = fsolve(eqc, [transformed_q0; transformed_psi0], options); % options defined above
    if flag <= 0
        transformed_q0_wtd = square_root_interval_to_real_line((1 - s.qmax_wt) * q0 + s.qmax * s.qmax_wt, s.qminA, s.qmax);
        [out1, ~, flag, ~] = fsolve(eqc, [transformed_q0_wtd; transformed_psi0], options);
        if flag <= 0
            opt2 = optimset('display', 'off', 'TolFun', 1e-6); % relax tolerance a bit
            [out1, ~, flag, ~] = fsolve(eqc, [transformed_q0; transformed_psi0], opt2);
            if flag <= 0
                % spline b/c q should be almost everywhere C2, but pchip for psi1 since not necessarily supposed to be C2
                q1   = interp1(grid.eta(1:j-1), grid.q(1:j-1), grid.eta(j), 'spline', 'extrap');
                q1   = min(q1, s.qmax);
                psi1 = interp1(grid.eta(1:j-1), grid.psiAa(1:j-1), grid.eta(j), 'pchip', 'extrap');
                out1 = [square_root_interval_to_real_line(q1, s.qminA, s.qmax); ...
                        square_root_interval_to_real_line(psi1, 0, 1)];
            end
        end
    end

    [~, grid] = eqcond(grid.eta(j), real(out1), fd_vals, v_interp, grid, s, 1);
    if grid.psiBa(j) <= 0
        break;
    end
    grid.qp(j) = (grid.q(j) - fd_vals(1)) / fd_vals(2);
    grid.flag(j) = flag;
    q0 = grid.q(j);
    psiAa0 = grid.psiAa(j);
    fd_vals = [q0; grid.eta(j+1) - grid.eta(j)];
end
new_start = j;

for j = new_start:n_eta
    eqc = @(x) eqcond(grid.eta(j), x, fd_vals, v_interp, grid, s, 2);

    % Transform initial guesses
    transformed_q0   = square_root_interval_to_real_line(q0, s.qminA, s.qmax);
    transformed_psi0 = square_root_interval_to_real_line(psiAa0, 0, 1);

    % start region 2, backward diff
    [out1, ~, flag, ~] = fsolve(eqc, [transformed_q0; transformed_psi0], options); % options defined above
    if flag <= 0
        transformed_q0_wtd = square_root_interval_to_real_line((1 - s.qmax_wt) * q0 + s.qmax * s.qmax_wt, s.qminA, s.qmax);
        [out1, ~, flag, ~] = fsolve(eqc, [transformed_q0_wtd; transformed_psi0], options);
        if flag <= 0
            opt2 = optimset('display', 'off', 'TolFun', 1e-6); % relax tolerance a bit
            [out1, ~, flag, ~] = fsolve(eqc, [transformed_q0; transformed_psi0], opt2);
            if flag <= 0
                % spline b/c q should be almost everywhere C2, but pchip for psi1 since not necessarily supposed to be C2
                q1   = interp1(grid.eta(floor(new_start/2):j-1), grid.q(floor(new_start/2):j-1), grid.eta(j), 'spline', 'extrap');
                q1   = min(q1, s.qmax);
                psi1 = interp1(grid.eta(floor(new_start/2):j-1), grid.psiAa(floor(new_start/2):j-1), grid.eta(j), 'pchip', 'extrap');
                out1 = [square_root_interval_to_real_line(q1, s.qminA, s.qmax); ...
                        square_root_interval_to_real_line(psi1, 0, 1)];
            end
        end
    end

    [~, grid] = eqcond(grid.eta(j), real(out1), fd_vals, v_interp, grid, s, 2);
    grid.qp(j) = (grid.q(j) - fd_vals(1)) / fd_vals(2);
    grid.flag(j) = flag;
    q0 = grid.q(j);
    psiAa0 = grid.psiAa(j);
    fd_vals = [q0; grid.eta(j+1) - grid.eta(j)];
end

%% Solve using forward differences from right side
s.side = 1;

% make guess for initial allocations and prices
eta0   = grid.eta(end);
psiBb0 = s.B_lvg * (1-eta0);
psiAb0 = (1 - (1 + (s.a / s.a_ * s.alpha / (1 - s.alpha))^s.s) * psiBb0) / ...
         (1 + (s.a / s.a_)^(s.s - 1) * (s.alpha / (1 - s.alpha))^s.s);
psiAa0 = 1 - psiBb0 - psiAb0;
Y = ces_production(s.a * psiAa0, s.a_ * psiAb0 + s.a * psiBb0, s.alpha, s.s);
q0 = (s.kappa * Y + 1) / (s.rho*s.kappa + s.lambda);
if s.psi ~= 1
    q0 = fsolve(@(q) get_q_eis(q, eta0, vA(end), vB(end), Y, s), q0, optimset('display', 'off'));
end

% Derive minimum capital price for right-hand side.
% The choice of psiAb is obtained from formula for psiAb
% when AE stops specializing by setting psiBb = 0
Ymin_psiAb = 1 / (1 + (s.a / s.a_)^(s.s - 1) * (s.alpha / (1 - s.alpha))^s.s);
Ymin_psiAa = 1 - Ymin_psiAb;
Ymin = ces_production(s.a * Ymin_psiAa, s.a_ * Ymin_psiAb, s.alpha, s.s);
if s.psi ~= 1
    qmin = fsolve(@(q) get_q_eis(q, 0, vA(end), vB(end), Ymin, s), q0, optimset('display', 'off'));
else
    qmin = (s.kappa * Ymin + 1) / (s.rho*s.kappa + s.lambda);
end
s.qminB = qmin; % if psi ~= 1, the qmin may be different when starting from B's side (left)

% start region 3, forward diff
fd_vals = [qmin; -(1 - grid.eta(end))];
for j = s.N:-1:(n_eta+1)
    eqc = @(x) eqcond(grid.eta(j), x, fd_vals, v_interp, grid, s, 3);

    % Transform initial guesses.
    transformed_q0   = square_root_interval_to_real_line(q0, s.qminB, s.qmax);
    transformed_psi0 = square_root_interval_to_real_line(psiBb0, 0, 1); % for side = 0, psiBb0 = x(2)

    % Solve!
    [out1, ~, flag, ~] = fsolve(eqc, [transformed_q0; transformed_psi0], options); % options defined above
    if flag <= 0
        transformed_q0_wtd = square_root_interval_to_real_line((1 - s.qmax_wt) * q0 + s.qmax * s.qmax_wt, s.qminB, s.qmax);
        [out1, ~, flag, ~] = fsolve(eqc, [transformed_q0_wtd; transformed_psi0], options);
        if flag <= 0
            opt2 = optimset('display', 'off', 'TolFun', 1e-6); % relax tolerance a bit
            [out1, ~, flag, ~] = fsolve(eqc, [transformed_q0; transformed_psi0], opt2);
            if flag <= 0
                % spline b/c q should be almost everywhere C2, but pchip for psi1 since not necessarily supposed to be C2
                q1   = interp1(grid.eta(j+1:s.N), grid.q(j+1:s.N), grid.eta(j), 'spline', 'extrap');
                q1   = min(q1, s.qmax);
                psi1 = interp1(grid.eta(j+1:s.N), grid.psiBb(j+1:s.N), grid.eta(j), 'pchip', 'extrap');
                out1 = [square_root_interval_to_real_line(q1, s.qminB, s.qmax); ...
                        square_root_interval_to_real_line(psi1, 0, 1)];
            end
        end
    end

    if numel(out1) > 2
        keyboard;
    end

    [~, grid] = eqcond(grid.eta(j), real(out1), fd_vals, v_interp, grid, s, 3);
    if grid.psiAb(j) <= 0
        break;
    end
    grid.qp(j) = (grid.q(j) - fd_vals(1)) / fd_vals(2);
    grid.flag(j) = flag;
    q0 = grid.q(j);
    psiBb0 = grid.psiBb(j);
    fd_vals = [q0; -(grid.eta(j) - grid.eta(j-1))];
end
new_start = j;

% start region 2, forward diff
for j = new_start:-1:(n_eta+1)
    eqc = @(x) eqcond(grid.eta(j), x, fd_vals, v_interp, grid, s, 2);

    % Transform initial guesses
    transformed_q0   = square_root_interval_to_real_line(q0, s.qminB, s.qmax);
    transformed_psi0 = square_root_interval_to_real_line(psiBb0, 0, 1);

    % start region 2, backward diff
    [out1, ~, flag, ~] = fsolve(eqc, [transformed_q0; transformed_psi0], options); % options defined above
    if flag <= 0
        transformed_q0_wtd = square_root_interval_to_real_line((1 - s.qmax_wt) * q0 + s.qmax * s.qmax_wt, s.qminB, s.qmax);
        [out1, ~, flag, ~] = fsolve(eqc, [transformed_q0_wtd; transformed_psi0], options);
        if flag <= 0
            opt2 = optimset('display', 'off', 'TolFun', 1e-6); % relax tolerance a bit
            [out1, ~, flag, ~] = fsolve(eqc, [transformed_q0; transformed_psi0], opt2);
            if flag <= 0
                % spline b/c q should be almost everywhere C2, but pchip for psi1 since not necessarily supposed to be C2
                q1   = interp1(grid.eta(j+1:new_start), grid.q(j+1:new_start), grid.eta(j), 'spline', 'extrap');
                q1   = min(q1, s.qmax);
                psi1 = interp1(grid.eta(j+1:new_start), grid.psiBb(j+1:new_start), grid.eta(j), 'pchip', 'extrap');
                out1 = [square_root_interval_to_real_line(q1, s.qminB, s.qmax); ...
                        square_root_interval_to_real_line(psi1, 0, 1)];
            end
        end
    end

    [~, grid] = eqcond(grid.eta(j), real(out1), fd_vals, v_interp, grid, s, 2);
    grid.qp(j) = (grid.q(j) - fd_vals(1)) / fd_vals(2);
    grid.flag(j) = flag;
    q0 = grid.q(j);
    psiBb0 = grid.psiBb(j);
    fd_vals = [q0; -(grid.eta(j) - grid.eta(j-1))];
end

% Compute central difference derivatives scheme
diff_coefs = get_fd_coef(2, 2, grid.eta, grid.eta(2:end-1));
diag_mat1 = [diff_coefs(2:end,1); 0; 0; 0];
diag_mat2 = [0; diff_coefs(2:end,2); 0; 0];
diag_mat3 = [0; 0; diff_coefs(2:end,3); 0];
diag_mat4 = [0; 0; 0; diff_coefs(2:end,4)];
diag_mat = [diag_mat1 diag_mat2 diag_mat3 diag_mat4];
D_xx = spdiags(diag_mat, -2:1, grid.dim1, grid.dim1);
D_xx(2, 1:4) = diff_coefs(1,:);

% compute q', q'', and other grid values
grid.qpp = D_xx * grid.q;
% grid = compute_vals(grid.eta, grid, vderivs, s);

clear v_interp; % Don't need to track these values anymore so clear from memory

%% Fill Grid values
if lastiter == 1
    % Compute risk free rate by interpolation
    qpp_interp        = griddedInterpolant(grid.eta(2:end-1), grid.qpp(2:end-1), 'pchip', 'extrap');
    grid.qpp(1)       = qpp_interp(grid.eta(1));
    grid.qpp(end)     = qpp_interp(grid.eta(end));
    grid.mu_q         = grid.qp ./ grid.q .* grid.mu_eta .* grid.eta + ...
                        grid.qpp ./ grid.q ./ 2  .* (grid.sig_etaA.^2 + ...
                        grid.sig_etaB.^2) .* grid.eta.^2;
    grid.iota         = 1 / s.kappa * (s.lambda * grid.q - 1);
    grid.rF           =  ((s.a .* grid.Pa - grid.iota) ./ ...
                        grid.q + grid.muK + grid.mu_q + ...
                        grid.sig_qA * s.sigA) - ...
						(s.gammaA .* grid.psiA ./ grid.eta .* ((s.sigA + ...
                        grid.sig_qA).^2 + grid.sig_qB.^2) + (s.gammaA-1) .* ...
                        ((s.sigA + grid.sig_qA) .* grid.sig_xiA + ...
                        grid.sig_qB .* grid.sig_xiB));
end

end

function out = find_argqmax(eta, v_interp, s)
    if isfield(s, 's') || s.s == 1
        psiA = s.alpha;
        psiB = 1 - s.alpha;
    else
        term1 = s.alpha^s.s;
        term2 = (1 - s.alpha)^s.s;
        psiA = term1 / (term1 + term2);
        psiB = term2 / (term1 + term2);
    end

    % finds the eta which maximizes q
    sigetaA = psiA * (1 / eta - 1) * s.sigA;
    sigetaB = -psiB * s.sigB;

    % central difference approximation at qmax
	vAp = v_interp.vAp_c(eta);
	vBp = v_interp.vBp_c(eta);
	vA  = v_interp.vA(eta);
	vB  = v_interp.vB(eta);
    sigvAA = vAp / vA * eta * sigetaA;
    sigvBB = vBp / vB * eta * sigetaB;
    sigxiA = sigvAA - sigetaA;
    sigzetaB = sigvBB + eta * sigetaB / (1 - eta);
    out = s.gammaA * psiA / eta * s.sigA^2 + ...
        (s.gammaA-1) * s.sigA * sigxiA - ...
        s.gammaB * psiB/ (1-eta) * s.sigB^2 - ...
        (s.gammaB-1) * s.sigB * sigzetaB;
end

function F = get_q_eis(q, eta, vA, vB, Y, s)
    if eta ~= 0 && eta ~= 1
        xi   = vA / eta / q;
        zeta = vB / (1 - eta) / q;
    elseif eta == 0
        zeta = vB / q;
        xi   = 1;
    elseif eta == 1
        xi   = vA / q;
        zeta = 1;
    end
    LHS  = s.kappa * (xi^(1-s.psi) * eta^s.psi + zeta^(1-s.psi) * (1-eta)^s.psi) * q^s.psi + s.lambda * q;
    F = LHS - s.kappa * Y - 1;
end

function [F, grid] = eqcond(eta, x, fd_vals, v_interp, grid, s, region)
% This function holds the implicit ODE that characterizes
% equilibrium in terms of the capital price q and its first
% derivative qp. This function also computes dynamics
%
% Written by William Chen, Gregory Phelan,  May 2021

% First element of x is always q
% and we transform it to always lie within [qmin, qmax]
if s.side == 0
    q = square_root_real_line_to_interval(x(1), s.qminA, s.qmax);
else
    q = square_root_real_line_to_interval(x(1), s.qminB, s.qmax);
end

%% Compute finite difference
qp = (q - fd_vals(1)) / fd_vals(2); % fd_vals(1) is q_{0} if forward diff, q_{-1} otherwise, fd_vals(2) is eta diff

%% Compute capital allocations

% Smoothly interpolate to get vA, vB
vA   = v_interp.vA(eta);
vB   = v_interp.vB(eta);
xi   = vA / eta / q;
% zeta = vB / (1 - eta) / q;

% Find PsiAs, PsiBs, first assuming specialization.
% x(2) is the guessed capital allocation variable,
% but it can be any one of psiAa, psiAb, etc.,
% depending on region
psi_val = square_root_real_line_to_interval(x(2), 0, 1);
lvgcons = 'NA'; % string to specify which agent is binding (NA = neither)
prat = s.a / s.a_; % productivity ratio
if region == 2
    % specialization
    psiAb = 0;
    psiBa = 0;
    if s.side == 0
        psiAa = psi_val;
    else
        psiAa = 1 - psi_val;
    end

    if psiAa / eta - 1 > s.LA
        lvgcons = 'A';
        psiAa = (1 + s.LA) * eta;
        psiBb = 1 - psiAa;
    elseif (1 - psiAa) / (1 - eta) - 1 > s.LB
        lvgcons = 'B';
        psiBb = (1 + s.LB) * (1 - eta);
        psiAa = 1 - psiBb;
    else
        psiBb = 1 - psiAa;
    end

    % Check specialization should hold. If not, change region
    if psiAa < s.threshold_cap_Aa
        region = 1;
    elseif psiBb < s.threshold_cap_Bb
        region = 3;
    end
end
if region == 1
    % B doesn't specialize
    psiAb = 0;
    psiAa = psi_val;
    psiBa = (1 - (1 + (prat * (1 - s.alpha) / s.alpha)^s.s) * psiAa) / ...
            (1 + (prat)^(s.s - 1) * ((1 - s.alpha) / s.alpha)^s.s);
    if psiAa / eta - 1 > s.LA
        lvgcons = 'A';
        psiAa = (1 + s.LA) * eta;
        psiBa = (1 - (1 + prat) * (1 + s.LA) * eta) / 2;
    end
    psiBb = 1 - psiAa - psiBa;
elseif region == 3
    % A doesn't specialize
    psiBa = 0;
    psiBb = psi_val;
    psiAb = (1 - (1 + (prat * s.alpha / (1 - s.alpha))^s.s) * psiBb) / ...
            (1 + (prat)^(s.s - 1) * (s.alpha / (1 - s.alpha))^s.s);
    if psiBb / (1 - eta) - 1 > s.LB
        lvgcons = 'B';
        psiBb = (1 - eta) * (1 + s.LB);
        psiAb = (1 - (1 + prat) * (1 - eta) * (1 + s.LB)) / 2;
    end
    psiAa = 1 - psiBb - psiAb;
end
ya   = s.a * psiAa + s.a_ * psiBa;
yb   = s.a * psiBb + s.a_ * psiAb;

% Compute output here. Needed to pin down capital allocations
% using the market-clearing condition for consumption
Y = ces_production(ya, yb, s.alpha, s.s);
%     if s.psi == 1
%         q = (s.kappa * output + 1) / (s.rho * s.kappa + s.lambda);
%     else
%         q = fsolve(@(qg) get_q_eis(qg, eta, vA, vB, output, s), q, optimset('display', 'off'));
%     end

% Compute prices and some allocations
Pa   = s.alpha * (Y / ya)^(1 / s.s); Pb = (1 - s.alpha) * (Y / yb)^(1 / s.s);
psiA = psiAa + psiAb;
psiB = psiBb + psiBa;

%% Compute volatilities
amplification = q - qp * (psiA - eta);

% Compute capital price volatilities
sig_qA = qp * psiA * (1 - eta) * s.sigA / amplification;
sig_qB = -qp * psiB * eta * s.sigB / amplification;
% sig_q = sqrt(sig_qA^2 + sig_qB^2);

% Compute eta volatilities
sig_etaA = psiA * (1 / eta - 1) * s.sigA + (psiA / eta - 1) * sig_qA;
sig_etaB = - psiB * s.sigB + (psiA / eta - 1) * sig_qB;
% sig_eta = sqrt(sigma_etaA^2 + sigma_etaB^2);

%% Loop to decide upwind direction
if region == 1 || region == 2
    loop_order = [1; 2; 3]; % guess forward difference first b/c expect mu_eta > 0
else
    loop_order = [2; 1; 3]; % guess backward difference first b/c expect mu_eta < 0
end
for loop_order_i = 1:3
    i = loop_order(loop_order_i);
    % Smoothly interpolate to get vAp
    if i == 1
        % First try forward difference
        vAp = v_interp.vAp_f(eta);
        vBp = v_interp.vBp_f(eta);
    elseif i == 2
        % Now backward difference
        vAp = v_interp.vAp_b(eta);
        vBp = v_interp.vBp_b(eta);
    else
        % Now central difference
        vAp = v_interp.vAp_c(eta);
        vBp = v_interp.vBp_c(eta);
    end

    % Compute volatilities of value functions
    sig_vAA = vAp / vA * sig_etaA * eta;
    sig_vAB = vAp / vA * sig_etaB * eta;
    sig_vBA = vBp / vB * sig_etaA * eta;
    sig_vBB = vBp / vB * sig_etaB * eta;

    % Compute volatilities of marginal value of wealth
    sig_xiA   = sig_vAA - (sig_etaA + sig_qA);
    sig_xiB   = sig_vAB - (sig_etaB + sig_qB);
    sig_zetaA = sig_vBA + eta / (1-eta) * sig_etaA - sig_qA;
    sig_zetaB = sig_vBB + eta / (1-eta) * sig_etaB - sig_qB;

    % Compute risk premia
    varsig_Aa    = psiA / eta * ((s.sigA + sig_qA)^2 + sig_qB^2);
    varsig_xiA   = (s.sigA + sig_qA) * sig_xiA + sig_qB * sig_xiB;
    varsig_Bb    = psiB / (1-eta) * (sig_qA^2 + (s.sigB + sig_qB)^2);
    varsig_zetaB = sig_qA * sig_zetaA + (s.sigB + sig_qB) * sig_zetaB;

    % Compute drift of eta
    vec_q  = [sig_qA; sig_qB];
    vec_K  = [psiA * s.sigA; psiB * s.sigB];
    vec_Aa = [s.sigA + sig_qA; sig_qB];
    % vec_Bb = [sig_qA; s.sigB + sig_qB]; % unnecessary to calculate this
    iota = 1 / s.kappa * (s.lambda * q - 1);
    if s.psi == 1
        c_nA = s.rho;
    else
        c_nA = xi^(1 - s.psi);
    end
    if isfield(s, 'eta_tau')
        indicator = eta < s.eta_tau;
    else
        indicator = 0;
    end
    if strcmp(lvgcons, 'NA')
        mu_eta = (s.a * Pa - iota) / q - c_nA + (psiA / eta - 1) * (s.gammaA * varsig_Aa + ...
            (s.gammaA - 1) * varsig_xiA) + (vec_q + vec_K)' * (vec_q + vec_K) + ...
            s.sigA * sig_qA - vec_q' * vec_K - (psiA / eta * vec_Aa)' * (vec_q + vec_K) + ...
            s.tau * (1 - eta) / eta * indicator;
    else
        mu_eta = (psiA / eta - 1) * (s.gammaB * varsig_Bb + (s.gammaB - 1) * varsig_zetaB - ...
            s.a * Pb / q) - iota / q + 1 / eta / q * (psiAa * s.a * Pa + ...
            psiAb * s.a_ * Pb) + (psiA / eta - 1 - psiA) * (s.sigA * sig_qA - ...
            s.sigB * sig_qB) - c_nA + (vec_q + vec_K)' * (vec_q + vec_K) - ...
            (psiA / eta * vec_Aa)' * (vec_q + vec_K) + s.tau * (1 - eta) / eta  * indicator;
    end

    if mu_eta > 0 && i == 1
        % Then forward difference is fine
        break;
    elseif mu_eta < 0 && i == 2
        % Then backward difference is fine
        break;
    end
end

%% Finish computations
% compute equilibrium asset pricing function
if strcmp(lvgcons, 'NA')
    % then use asset pricing condition to pin down q
	% and market-clearing for capital to pin down psi_val
    returns_diff   = s.a * (Pa - Pb) / q + s.sigA * sig_qA - s.sigB * sig_qB;
    risk_prem_diff = s.gammaA * varsig_Aa + (s.gammaA - 1) * varsig_xiA - ...
                     s.gammaB * varsig_Bb - (s.gammaB - 1) * varsig_zetaB;
    mc_consumption = Y - s.rho * q - iota;

    F              = [returns_diff - risk_prem_diff; mc_consumption];
else
    % then use market clearing for consumption to pin down q
    % and use the leverage constraint equation for psi_val.
    % Note that we do not use Y - s.rho * q - iota
    % b/c iota is computed based on the guessed q value.
    % What we want to do for the error is compute the
    % implied q based on the fact that we know capital
    % allocations given binding leverage constraints
    % and then compare the implied q to the guessed on.
    if s.psi == 1
        q_implied = (s.kappa * Y + 1) / (s.rho * s.kappa + s.lambda);
    else
        q_implied = fsolve(@(qg) get_q_eis(qg, eta, vA, vB, Y, s), q, optimset('display', 'off'));
    end
    if strcmp(lvgcons, 'A')
        psi_implied = psiAa; % = (1 + s.LA) * eta;
    else
        psi_implied = psiBb; % = (1 + s.LB) * (1 - eta);
    end

    F = [q - q_implied; psi_val - psi_implied];
end

%% Compute grid values
if nargout > 1
    i = find(grid.eta == eta);
    muK   = s.lambda / s.kappa * log(s.lambda * q) - s.delta;

    % Store variables
    grid.q(i)         = q;
    grid.mu_eta(i)    = mu_eta;
    grid.sig_etaA(i)  = sig_etaA;
    grid.sig_etaB(i)  = sig_etaB;
    grid.psiAa(i)     = real(psiAa);
    grid.psiAb(i)     = real(psiAb);
    grid.psiBa(i)     = real(psiBa);
    grid.psiBb(i)     = real(psiBb);
    grid.psiA(i)      = real(psiA);
    grid.psiB(i)      = real(psiB);
    grid.sig_vAA(i)   = sig_vAA;
    grid.sig_vAB(i)   = sig_vAB;
    grid.sig_vBA(i)   = sig_vBA;
    grid.sig_vBB(i)   = sig_vBB;
    grid.muK(i)       = muK;
    grid.sig_qA(i)    = sig_qA;
    grid.sig_qB(i)    = sig_qB;
    grid.sig_xiA(i)   = sig_xiA;
    grid.sig_xiB(i)   = sig_xiB;
    grid.Pa(i)        = Pa;
    grid.Pb(i)        = Pb;
    grid.sig_zetaA(i) = sig_zetaA;
    grid.sig_zetaB(i) = sig_zetaB;
end
end

function out = ces_production(ya, yb, alpha, s)
% Computes the CES output given intermediates ya, yb
% with weighting alpha and CES s
    if s == 1
        out = ya^alpha * yb^(1-alpha);
    else
        out_pwr = s / (s - 1); % power on the outside
        in_pwr  = 1 / out_pwr; % power on the inside
        out = (alpha * ya^in_pwr + (1 - alpha) * yb^in_pwr)^out_pwr;
    end
end
