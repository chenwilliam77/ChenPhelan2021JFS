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
qmax = (s.kappa * max(s.aA, s.aB) + 1) / (s.rho * s.kappa + s.lambda);
s.qmax = qmax;
options = optimset('display', 'off', 'TolFun', s.static_tol);
initpt = s.alpha * (s.gammaB / (s.gammaA + s.gammaB));
eta_maxq = fsolve(@(eta) find_argqmax(eta, grid, vA, vB, vAp_c(2:end-1), vBp_c(2:end-1), s), initpt, options);

%% Solve from the left first
% Make guess for capital allocation
eta0   = grid.eta(1);
psiA0 = s.A_lvg * eta0;
psiB0 = 1 - psiA0;
Y = psiA0 * s.aA + psiB0 * s.aB;
q0 = (s.kappa * Y + 1)/(s.rho * s.kappa + s.lambda);
if s.psi ~= 1
    q0 = fsolve(@(q) get_q_eis(q, eta0, vA(1), vB(1), Y, s), q0, optimset('display', 'off'));
    qmin = fsolve(@(q) get_q_eis(q, 0, vA(1), vB(1), min(s.aA, s.aB), s), q0, optimset('display', 'off'));
else
    qmin = (s.kappa * min(s.aA, s.aB) + 1)/(s.rho*s.kappa + s.lambda);
end
grid.qminA = qmin;

% eta0 = interp1(grid.eta, grid.eta, eta_maxq, 'nearest');
% n_eta = find(grid.eta == eta0); % location of the initial starting eta, which maximizes q
% if eta0 == eta_maxq
%     grid.eta(n_eta) = eta0 - (eta0 - grid.eta(n_eta-1))/3;
%     % eta0 = grid.eta(n_eta);
% elseif eta0 > eta_maxq
%     n_eta = n_eta - 1;
%     % eta0 = grid.eta(n_eta);
% end

%% Solve using backward differences
% start with region 1, backward diff
fd_vals = [qmin; s.start];
grid.flag = zeros(size(grid.eta));
for j = 1:grid.dim1
    eqc = @(q) eqcond(grid.eta(j), q, fd_vals, vderivs, grid, s, 1);
    [q1, ~, flag, ~] = fsolve(eqc, q0, options); % options defined above
    q1 = real(q1);
    [~, grid] = eqcond(grid.eta(j), q1, fd_vals, vderivs, grid, s, 1);
    grid.q(j) = q1;
    grid.qp(j) = (q1 - fd_vals(1)) / fd_vals(2);
    grid.flag(j) = flag;
    q0 = q1;
    if j < grid.dim1
        fd_vals = [q0; grid.eta(j+1) - grid.eta(j)];
    end
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
    grid.qpp(2:end-1) = (grid.q(3:end) - grid.q(1:end-2)) ./ (grid.eta(3:end) - grid.eta(1:end-2));
    grid.qpp(1)       = interp1(grid.eta(2:end-1), grid.qpp(2:end-1), grid.eta(end), 'pchip', 'extrap');
    grid.qpp(end)     = interp1(grid.eta(2:end-1), grid.qpp(2:end-1), grid.eta(end), 'pchip', 'extrap');
    grid.mu_q         = grid.qp ./ grid.q .* grid.mu_eta .* grid.eta + ...
                        grid.qpp ./ grid.q ./ 2  .* (grid.sig_etaA.^2 + ...
                        grid.sig_etaB.^2) .* grid.eta.^2;
    grid.iota         = 1 / s.kappa * (s.lambda * grid.q - 1);
    grid.rF           = ((s.aA - grid.iota) ./ grid.q + grid.muK + grid.mu_q + ...
                        grid.sig_qA * s.sigA) - s.gammaA .* grid.psiA ./ grid.eta .* ((s.sigA + ...
                        grid.sig_qA).^2 + grid.sig_qB.^2) + (s.gammaA-1) .* ...
                        ((s.sigA + grid.sig_qA) .* grid.sig_xiA + ...
                        grid.sig_qB .* grid.sig_xiB);
end

end

function out = find_argqmax(eta, grid, vAvec, vBvec, vApvec, vBpvec, s)
    % finds the eta which maximizes q
    sigetaA = 1 / 2 * (1 / eta - 1) * s.sigA;
    sigetaB = 1 / 2 * (1/eta-1) * s.sigB;

    % central difference approximation at qmax
    vAp = interp1(grid.eta(2:end-1), vApvec, eta, 'pchip', 'extrap');
    vBp = interp1(grid.eta(2:end-1), vBpvec, eta, 'pchip', 'extrap');
    vA = interp1(grid.eta, vAvec, eta, 'pchip', 'extrap');
    vB = interp1(grid.eta, vBvec, eta, 'pchip', 'extrap');
    sigvAA = vAp / vA * eta * sigetaA;
    sigvBB = vBp / vB * eta * sigetaB;
    sigxiA = sigvAA + sigetaA;
    sigzetaB = sigvBB + sigetaB;
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

% Find PsiAs, PsiBs, first assuming specialization.  This assumes Cobb-Douglas
% with even weights, otherwise you need to use fsolve.
lvgcons = 0; % indicator for whether an agent is leverage constrained
if region == 2
    % country A produces all
    psiA = 1;
    psiB = 1 - psiA;
elseif region == 1
    % B holds some capital
    psiA = ((s.rho * s.kappa + s.lambda) * q - 1 - s.kappa * s.aB) / (s.aA - s.aB) / s.kappa;
    if psiA / eta - 1 > s.LA
        lvgcons = 1;
        psiA = (1 + s.LA) * eta;
    elseif (1 - psiA) / (1 - eta) - 1 > s.LB
        lvgcons = 1;
        psiA = 1 - (1 + s.LB) * (1 - eta);
    end
    psiB = 1 - psiA;
end
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
    if lvgcons == 0
        mu_eta = (s.aA - iota) / q - c_nA + (psiA / eta - 1) * (s.gammaA * varsig_Aa + ...
            (s.gammaA - 1) * varsig_xiA) + (vec_q + vec_K)' * (vec_q + vec_K) + ...
            s.sigA * sig_qA - vec_q' * vec_K - (psiA / eta * vec_Aa)' * (vec_q + vec_K) + s.tau * (1 - eta) / eta - s.nu;
    else
        mu_eta = (psiA / eta - 1) * (s.gammaB * varsig_Bb + (s.gammaB - 1) * varsig_zetaB - ...
            s.aB / q) - iota / q + 1 / eta / q * (psiA * s.aA) + ...
            (psiA / eta - 1 - psiA) * (s.sigA * sig_qA - ...
            s.sigB * sig_qB) - c_nA + (vec_q + vec_K)' * (vec_q + vec_K) - ...
            (psiA / eta * vec_Aa)' * (vec_q + vec_K) + s.tau * (1 - eta) / eta - s.nu;
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
returns_diff   = (s.aA - s.aB) / q + s.sigA * sig_qA - s.sigB * sig_qB;
risk_prem_diff = s.gammaA * varsig_Aa + (s.gammaA - 1) * varsig_xiA - ...
                    s.gammaB * varsig_Bb - (s.gammaB - 1) * varsig_zetaB;
if lvgcons == 0
    % then use asset pricing difference to pin down q
    F = returns_diff - risk_prem_diff;
elseif psiA >= 1
    F = q - s.qmax;
else
    % then use market clearing for consumption to pin down q
    Y = s.aA * psiA + s.aB * psiB;
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
    grid.sig_zetaA(i) = sig_zetaA;
    grid.sig_zetaB(i) = sig_zetaB;
end
end
