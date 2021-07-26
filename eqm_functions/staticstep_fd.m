function [grid,s] = staticstep_fd(grid, vA, vB, s, lastiter)
% This function performs the static step to solve for equilibrium,
% given vA and vB
%
% Written by William Chen, Mar. 2019

if nargin < 5
    lastiter = 0;
elseif not(lastiter == 1 || lastiter == 0)
    error('Last argument is not 0 or 1');
end

%% Set up
% Compute first-order derivatives of vA, vB for equilibrium computations
deta = diff(grid.eta); vAp = diff(vA) ./ deta; vBp = diff(vB) ./ deta;
vAp_f = [vAp; 0]; vAp_b = [0; vAp]; vBp_f = [vBp; 0]; vBp_b = [0; vBp];
vAp_c = [0; (vA(3:end) - vA(1:end-2)) ./ sum(deta(2:end) + deta(1:end-1)); 0];
vBp_c = [0; (vB(3:end) - vB(1:end-2)) ./ sum(deta(2:end) + deta(1:end-1)); 0];

if ~isreal(vA) || ~isreal(vB)
    keyboard;
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

% Store derivatives
vderivs.eta = grid.eta; vderivs.vA = vA; vderivs.vB = vB;
vderivs.vAp_f = vAp_f; vderivs.vBp_f = vBp_f;
vderivs.vAp_b = vAp_b; vderivs.vBp_b = vBp_b;
vderivs.vAp_c = vAp_c; vderivs.vBp_c = vBp_c;
vderivs.vApp      = D_xx * vA;
vderivs.vBpp      = D_xx * vB;
vderivs.vApp(1)   = interp1(grid.eta(2:end-1), vderivs.vApp(2:end-1), grid.eta(1), 'pchip', 'extrap');
vderivs.vApp(end) = interp1(grid.eta(2:end-1), vderivs.vApp(2:end-1), grid.eta(end), 'pchip', 'extrap');
vderivs.vBpp(1)   = interp1(grid.eta(2:end-1), vderivs.vBpp(2:end-1), grid.eta(1), 'pchip', 'extrap');
vderivs.vBpp(end) = interp1(grid.eta(2:end-1), vderivs.vBpp(2:end-1), grid.eta(end), 'pchip', 'extrap');

% Find eta which maximizes capital price
qmax = (s.kappa * s.a * (s.alpha)^s.alpha * (1-s.alpha)^(1-s.alpha) + 1) / ...
        (s.rho * s.kappa + s.lambda);
s.qmax = qmax;
options = optimset('display', 'off', 'TolFun', s.static_tol);
initpt = s.alpha * (s.gammaB / (s.gammaA + s.gammaB));

% Find minimum capital allocation before specialization stops
s.max_cap = 1 / (1 + s.a / s.a_ * (1 - s.alpha) / s.alpha);


%% Solve from the left first
% Make guess for capital allocation
eta0   = grid.eta(1);
psiAa0 = s.A_lvg * eta0;
psiBa0 = (1- (1 + s.a / s.a_) * psiAa0) / 2;
psiBb0 = s.a / s.a_ * psiAa0 + psiBa0;
Y = sqrt((s.a * psiAa0 + s.a_ * psiBa0) * s.a * psiBb0);
q0 = (s.kappa * Y + 1)/(s.rho * s.kappa + s.lambda);
if s.psi ~= 1
    q0 = fsolve(@(q) get_q_eis(q, eta0, vA(1), vB(1), Y, s), q0, optimset('display', 'off'));
    qmin = fsolve(@(q) get_q_eis(q, 0, vA(1), vB(1), sqrt(s.a * s.a_) / 2, s), q0, optimset('display', 'off'));
else
    qmin = (s.kappa * (sqrt(s.a_ * s.a) / 2) + 1)/(s.rho*s.kappa + s.lambda);
end
grid.qminA = qmin;

% n_eta = floor(s.N / 2);
eta_maxq = fsolve(@(eta) find_argqmax(eta, grid, vA, vB, vAp_c(2:end-1), vBp_c(2:end-1), s), initpt, options);
eta0 = interp1(grid.eta, grid.eta, eta_maxq, 'nearest');
n_eta = find(grid.eta == eta0); % location of the initial starting eta, which maximizes q
if eta0 == eta_maxq
    grid.eta(n_eta) = eta0 - (eta0 - grid.eta(n_eta-1))/3;
    % eta0 = grid.eta(n_eta);
elseif eta0 > eta_maxq
    n_eta = n_eta - 1;
    % eta0 = grid.eta(n_eta);
end
if s.gammaA == s.gammaB
    n_eta = ceil(s.N / 2);
else
    n_eta = ceil(2/3 * n_eta + 1/3 * s.N / 2);
end



%% Solve using backward differences
% start with region 1, backward diff
fd_vals = [qmin; s.start];
grid.flag = zeros(size(grid.eta));
s.side = 0; % starting before max q
for j = 1:n_eta
    eqc = @(q) eqcond(grid.eta(j), q, fd_vals, vderivs, grid, s, 1);
    [q1, ~, flag, ~] = fsolve(eqc, q0, options); % options defined above
    if flag <= 0
        [q1, ~, flag, ~] = fsolve(eqc, q0 * (1 - s.qmax_wt) + qmax * s.qmax_wt, options);
        if flag <= 0
            opt2 = optimset('display', 'off', 'TolFun', 1e-6);
            [q1, ~, flag, ~] = fsolve(eqc, q0 * (1 - s.qmax_wt) + qmax * s.qmax_wt, opt2);
            if flag <= 0
                q1 = interp1(grid.eta(1:j-1), grid.q(1:j-1), grid.eta(j), 'spline', 'extrap');
            end
        end
    end
    q1 = real(q1);
    [~, grid] = eqcond(grid.eta(j), q1, fd_vals, vderivs, grid, s, 1);
    if grid.psiBa(j) <= 0
        break;
    end
    grid.q(j) = q1;
    grid.qp(j) = (q1 - fd_vals(1)) / fd_vals(2);
    grid.flag(j) = flag;
    q0 = q1;
    fd_vals = [q0; grid.eta(j+1) - grid.eta(j)];
end

if j < n_eta
    new_start = j;
    % start region 2, backward diff
    for j = new_start:n_eta
        eqc = @(q) eqcond(grid.eta(j), q, fd_vals, vderivs, grid, s, 2);
        [q1, ~, flag, ~] = fsolve(eqc, q0, options);
        if flag <= 0
            [q1, ~, flag, ~] = fsolve(eqc, q0 * (1 - s.qmax_wt) + qmax * s.qmax_wt, options);
            if flag <= 0
                opt2 = optimset('display', 'off', 'TolFun', 1e-6);
                [q1, ~, flag, ~] = fsolve(eqc, q0 * (1 - s.qmax_wt) + qmax * s.qmax_wt, opt2);
                if flag <= 0
                    q1 = interp1(grid.eta(floor(new_start/2):j-1), grid.q(floor(new_start/2):j-1), grid.eta(j), 'spline', 'extrap');
                end
            end
        end
        q1 = real(q1);
        [~, grid] = eqcond(grid.eta(j), q1, fd_vals, vderivs, grid, s, 2);
        grid.q(j) = q1;
        grid.qp(j) = (q1 - fd_vals(1)) / fd_vals(2);
        grid.flag(j) = flag;
        q0 = q1;
        fd_vals = [q0; grid.eta(j+1) - grid.eta(j)];
    end
end

%% Solve using forward differences from right side
s.side = 1;

% make guess
eta0   = grid.eta(end);
psiBb0 = s.B_lvg * (1-eta0);
psiAb0 = (1- (1 + s.a / s.a_) * psiBb0) / 2;
psiAa0 = s.a / s.a_ * psiBb0 + psiAb0;
Y  = sqrt((s.a * psiAa0) * (s.a * psiBb0 + s.a_ * psiAb0));
q0 = (s.kappa * Y + 1)/(s.rho*s.kappa + s.lambda);
if s.psi ~= 1
    q0 = fsolve(@(q) get_q_eis(q, eta0, vA(end), vB(end), Y, s), q0, optimset('display', 'off'));
    qmin = fsolve(@(q) get_q_eis(q, 1, vA(end), vB(end), sqrt(s.a * s.a_) / 2, s), q0, optimset('display', 'off'));
else
    qmin = (s.kappa * (sqrt(s.a_ * s.a) / 2) + 1)/(s.rho * s.kappa + s.lambda);
end
grid.qminB = qmin;

% start region 3, forward diff
fd_vals = [qmin; -(1 - grid.eta(end))];
for j = length(grid.eta):-1:(n_eta+1)
    eqc = @(q) eqcond(grid.eta(j), q, fd_vals, vderivs, grid, s, 3);
    [q1, ~, flag, ~] = fsolve(eqc, q0 * (1 - s.qmax_wt) + qmax * s.qmax_wt, options);
    if flag <= 0
        [q1, ~, flag, ~] = fsolve(eqc, q0, options);
        if flag <= 0
            opt2 = optimset('display', 'off', 'TolFun', 1e-6);
            [q1, ~, flag, ~] = fsolve(eqc, q0 * (1 - s.qmax_wt) + qmax * s.qmax_wt, opt2);
            if flag <= 0
                q1 = interp1(grid.eta(j+1:length(grid.eta)), grid.q(j+1:length(grid.eta)), grid.eta(j), 'spline', 'extrap');
            end
        end
    end
    if numel(q1) > 1
        keyboard;
    end
    q1 = real(q1);
    [~, grid] = eqcond(grid.eta(j), q1, fd_vals, vderivs, grid, s, 3);
    if grid.psiAb(j) < 0
        break;
    end
    grid.q(j) = q1;
    grid.qp(j) = (q1 - fd_vals(1)) / fd_vals(2);
    grid.flag(j) = flag;
    q0 = q1;
    fd_vals = [q0; -(grid.eta(j) - grid.eta(j-1))];
end
new_start = j;

% start region 2, forward diff
for j = new_start:-1:(n_eta+1)
    eqc = @(q) eqcond(grid.eta(j), q, fd_vals, vderivs, grid, s, 2);
    [q1, ~, flag, ~] = fsolve(eqc, q0, options);
    if flag <= 0
        [q1, ~, flag, ~] = fsolve(eqc, q0 * (1 - s.qmax_wt) + qmax * s.qmax_wt, options);
        if flag <= 0
            opt2 = optimset('display', 'off', 'TolFun', 1e-6);
            [q1, ~, flag, ~] = fsolve(eqc, q0 * (1 - s.qmax_wt) + qmax * s.qmax_wt, opt2);
            if flag <= 0
                up_ind = floor((length(grid.eta) - (n_eta + 1)) / 2);
                q1 = interp1(grid.eta(j+1:j+up_ind), grid.q(j+1:j+up_ind), grid.eta(j), 'spline', 'extrap');
                q1 = min(q1, qmax);
            end
        end
    end
    q1 = real(q1);
    [~, grid] = eqcond(grid.eta(j), q1, fd_vals, vderivs, grid, s, 2);
    grid.q(j) = q1;
    grid.qp(j) = (q1 - fd_vals(1)) / fd_vals(2);
    grid.flag(j) = flag;
    q0 = q1;
    fd_vals = [q0; -(grid.eta(j) - grid.eta(j-1))];
end

% compute q', q'', and other grid values
grid.qpp = D_xx * grid.q;
% grid = compute_vals(grid.eta, grid, vderivs, s);

clear vderivs; % Don't need to track these values anymore so clear from memory

%% Fill Grid values
% grid.q = grid.q;
% grid.qp = grid.qp;
% grid.mu_eta   = grid.mu_eta;
% grid.sig_etaA = grid.sig_etaA;
% grid.sig_etaB = grid.sig_etaB;
% grid.PsiA     = grid.PsiA;
% grid.PsiB     = grid.PsiB;
% grid.sig_vAA  = grid.sig_vAA;
% grid.sig_vAB  = grid.sig_vAB;
% grid.sig_vBA  = grid.sig_vBA;
% grid.sig_vBB  = grid.sig_vBB;
% grid.muK      = grid.muK;
if lastiter == 1
%     grid.sig_qA    = grid.sig_qA;
%     grid.sig_qB    = grid.sig_qB;
%     grid.sig_xiA   = grid.sig_xiA;
%     grid.sig_xiB   = grid.sig_xiB;
%     grid.sig_zetaA = grid.sig_zetaA;
%     grid.sig_zetaB = grid.sig_zetaB;
%     grid.PsiAa     = grid.PsiAa;
%     grid.PsiBa     = grid.PsiBa;
%     grid.PsiAb     = grid.PsiAb;
%     grid.PsiBb     = grid.PsiBb;
%     grid.Pa        = grid.Pa;
%     grid.Pb        = grid.Pb;

    % Compute risk free rate by interpolation
    grid.qpp(1)       = interp1(grid.eta(2:end-1), grid.qpp(2:end-1), grid.eta(end), 'pchip', 'extrap');
    grid.qpp(end)     = interp1(grid.eta(2:end-1), grid.qpp(2:end-1), grid.eta(end), 'pchip', 'extrap');
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

function out = find_argqmax(eta, grid, vAvec, vBvec, vApvec, vBpvec, s)
    % finds the eta which maximizes q
    sigetaA = 1 / 2 * (1 / eta - 1) * s.sigA;
    sigetaB = 1 / 2 * (1/eta-1) * s.sigB;
    % The line above for sigetaB is technically a bug.
    % It really should sigetaB = -(1 / 2) * s.sigB;
    % However, this function is just used to find a good eta
    % at which we want to split the state space into the two "sides.
    % Thus, this bug just means that we may be splitting the state space
    % into two regions at a different point than anticipated.
    % The numbers for things like ExpVA are slightly different
    % due to the implied change in the approximation scheme, but it
    % does not affect our main results e.g. Nash is LA = LB = 0;
    % the AE prefers Nash to CE, unlike the EME; etc.
    % This bug is left here to avoid having to regenerate all the results
    % with slightly different numbers.

    % central difference approximation at qmax
    vAp = interp1(grid.eta(2:end-1), vApvec, eta, 'pchip', 'extrap');
    vBp = interp1(grid.eta(2:end-1), vBpvec, eta, 'pchip', 'extrap');
    vA = interp1(grid.eta, vAvec, eta, 'pchip', 'extrap');
    vB = interp1(grid.eta, vBvec, eta, 'pchip', 'extrap');
    sigvAA = vAp / vA * eta * sigetaA;
    sigvBB = vBp / vB * eta * sigetaB;
    sigxiA = sigvAA - sigetaA;
    sigzetaB = sigvBB + eta * sigetaB / (1 - eta);
    out = s.gammaA * 1 / (2 * eta) * s.sigA^2 + ...
        (s.gammaA-1) * s.sigA * sigxiA - ...
        s.gammaB * 1 / (2 * (1-eta)) * s.sigB^2 - ...
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

function [F, grid] = eqcond(eta, q, fd_vals, vderivs, grid, s, region)
% This function holds the implicit ODE that characterizes
% equilibrium in terms of the capital price q and its first
% derivative qp. This function also computes dynamics
%
% Written by William Chen, Gregory Phelan Mar. 2019

if q > s.qmax
    F = 1e2;
    return
end
%% Compute finite difference
qp = (q - fd_vals(1)) / fd_vals(2); % fd_vals(1) is q_{0} if forward diff, q_{-1} otherwise, fd_vals(2) is eta diff

%% Compute capital allocations

% Smoothly interpolate to get vA, vB
vA   = interp1(vderivs.eta, vderivs.vA, eta, 'spline', 'extrap');
vB   = interp1(vderivs.eta, vderivs.vB, eta, 'spline', 'extrap');
xi   = vA / eta / q;
zeta = vB / (1 - eta) / q;

% Compute output
if s.psi == 1
    Y = ((s.rho * s.kappa + s.lambda) * q - 1) / s.kappa;
else
    Y = ((s.kappa * (xi^(1 - s.psi) * eta + zeta^(1 - s.psi) * (1 - eta)) + ...
         s.lambda) * q - 1) / s.kappa;
end

% Find PsiAs, PsiBs, first assuming specialization.  This assumes Cobb-Douglas
% with even weights, otherwise you need to use fsolve.
lvgcons = 0; % indicator for whether an agent is leverage constrained
if region == 2
    % specialization
    psiAb = 0;
    psiBa = 0;
    Cterm = Y / s.a;
    % psiAa = x(2);
    if s.side == 0
        psiAa = (1 - sqrt(1 - 4 * Cterm^2)) / 2;
    else
        psiAa = (1 + sqrt(1 - 4 * Cterm^2)) / 2;
    end
    if psiAa / eta - 1 > s.LA
        lvgcons = 1;
        psiAa = (1 + s.LA) * eta;
        psiBb = 1 - psiAa;
    elseif (1 - psiAa) / (1 - eta) - 1 > s.LB
        lvgcons = 1;
        psiBb = (1 + s.LB) * (1 - eta);
        psiAa = 1 - psiBb;
    else
        psiBb = 1 - psiAa;
    end

    % Check specialization should hold. If not, change region
    if psiAa < s.max_cap
        region = 1;
    elseif psiBb < s.max_cap
        region = 3;
    end
end
if region == 1
    % B doesn't specialize
    prat = s.a / s.a_;
    psiAb = 0;
    psiBa = (prat / (prat - 1)) * (1 - (1 + prat) * Y / (s.a * sqrt(prat)));
    psiAa = (1 - 2 * psiBa) / (1 + prat);
    if psiAa / eta - 1 > s.LA
        lvgcons = 1;
        psiAa = (1 + s.LA) * eta;
        psiBa = (1 - (1 + prat) * (1 + s.LA) * eta) / 2;
        psiBb = 1 - psiAa - psiBa;
%     elseif (1 - psiAa) / (1 - eta) - 1 > s.LB
%         lvgcons = 1;
%         psiAa = 1 - (1 + s.LB) * (1 - eta);
%         psiBa = (1 - (1 + prat) * psiAa) / 2;
%         psiBb = 1 - psiAa - psiBa;
    else
        psiBb = prat * psiAa + psiBa;
    end
elseif region == 3
    % A doesn't specialize
    prat = s.a / s.a_;
    psiBa = 0;
    psiAb = (prat / (prat - 1)) * (1 - (1 + prat) * Y / (s.a * sqrt(prat)));
    psiBb = (1 - 2 * psiAb) / (1 + prat);
    if psiBb / (1 - eta) - 1 > s.LB
        lvgcons = 1;
        psiBb = (1 - eta) * (1 + s.LB);
        psiAb = (1 - (1 + prat) * (1 - eta) * (1 + s.LB)) / 2;
        psiAa = 1 - psiBb - psiAb;
    else
        psiAa =  prat * psiBb + psiAb;
    end
end
ya   = s.a * psiAa + s.a_ * psiBa;
yb   = s.a * psiBb + s.a_ * psiAb;
if lvgcons == 1
    % Need to choose correct output
    Y = ya^s.alpha * yb^(1 - s.alpha);
%     if s.psi == 1
%         q = (s.kappa * output + 1) / (s.rho * s.kappa + s.lambda);
%     else
%         q = fsolve(@(qg) get_q_eis(qg, eta, vA, vB, output, s), q, optimset('display', 'off'));
%     end
end
Pa   = s.alpha * (Y/ya); Pb = (1 - s.alpha) * (Y/yb);
psiA = psiAa + psiAb;
psiB = psiBb + psiBa;
amplification = q - qp * (psiA - eta);

%% Compute volatilities
% Compute capital price volatilities
sig_qA = qp * psiA * (1 - eta) * s.sigA / amplification;
sig_qB = -qp * psiB * eta * s.sigB / amplification;
% sig_q = sqrt(sig_qA^2 + sig_qB^2);

% Compute eta volatilities
sig_etaA = psiA * (1 / eta - 1) * s.sigA + (psiA / eta - 1) * sig_qA;
sig_etaB = - psiB * s.sigB + (psiA / eta - 1) * sig_qB;
% sig_eta = sqrt(sigma_etaA^2 + sigma_etaB^2);

%% Loop to decide upwind direction
for i = 1:3
    % Smoothly interpolate to get vAp
    if i == 1
        % First try forward difference
        vAp = interp1(vderivs.eta(1:end-1), vderivs.vAp_f(1:end-1), eta, 'spline', 'extrap');
        vBp = interp1(vderivs.eta(1:end-1), vderivs.vBp_f(1:end-1), eta, 'spline', 'extrap');
    elseif i == 2
        % Now backward difference
        vAp = interp1(vderivs.eta(2:end), vderivs.vAp_b(2:end), eta, 'spline', 'extrap');
        vBp = interp1(vderivs.eta(2:end), vderivs.vBp_b(2:end), eta, 'spline', 'extrap');
    else
        % Now central difference
        vAp = interp1(vderivs.eta(2:end-1), vderivs.vAp_c(2:end-1), eta, 'spline', 'extrap');
        vBp = interp1(vderivs.eta(2:end-1), vderivs.vBp_c(2:end-1), eta, 'spline', 'extrap');
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
    % vec_Bb = [sig_qA; s.sigB + sig_qB];
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
    if lvgcons == 0
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
returns_diff   = s.a * (Pa - Pb) / q + s.sigA * sig_qA - s.sigB * sig_qB;
risk_prem_diff = s.gammaA * varsig_Aa + (s.gammaA - 1) * varsig_xiA - ...
                    s.gammaB * varsig_Bb - (s.gammaB - 1) * varsig_zetaB;
if lvgcons == 0
    % then use asset pricing difference to pin down q
    F = returns_diff - risk_prem_diff;
else
    % then use market clearing for consumption to pin down q
    if s.psi == 1
        qL = (s.kappa * Y + 1) / (s.rho * s.kappa + s.lambda);
    else
        qL = fsolve(@(qg) get_q_eis(qg, eta, vA, vB, Y, s), q, optimset('display', 'off'));
    end
    F = q - qL;
end

%% Compute grid values
if nargout > 1
    i = find(grid.eta == eta);
    muK   = s.lambda / s.kappa * log(s.lambda * q) - s.delta;
    vApp  = interp1(grid.eta, vderivs.vApp, eta, 'nearest');
    vBpp  = interp1(grid.eta, vderivs.vBpp, eta, 'nearest');
    mu_vA1 = s.gammaA / 2 * ((psiA * s.sigA + sig_vAA).^2 + ...
            (psiB * s.sigB + sig_vAB).^2) + s.rho * (log(vA) - ...
            log(s.rho * q .* eta)) - muK - ...
            psiA .* s.sigA .* sig_vAA - ...
            psiB * s.sigB * sig_vAB;
    mu_vA2 = vAp / vA * mu_eta * eta + vApp / vA * eta^2 * (sig_etaA^2 + sig_etaB^2);
    mu_vB1 = s.gammaB / 2 * ((psiA * s.sigA + sig_vBA).^2 + ...
            (psiB * s.sigB + sig_vBB).^2) + s.rho * (log(vB) - ...
            log(s.rho * q .* (1 - eta))) - ...
            muK - psiA * s.sigA * sig_vBA - ...
            psiB * s.sigB * sig_vBB;
    mu_vB2 = vBp / vB * mu_eta * eta + vBpp / vB * eta^2 * (sig_etaA^2 + sig_etaB^2);
    mu_vAerr = mu_vA1 - mu_vA2;
    mu_vBerr = mu_vB1 - mu_vB2;

    % Store variables
    grid.mu_eta(i)    = mu_eta;
    grid.mu_vAerr(i)  = mu_vAerr;
    grid.mu_vBerr(i)  = mu_vBerr;
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
