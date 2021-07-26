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

% Set capital price
grid.q = ones(size(grid.eta)) / s.lambda;
options = optimset('display', 'off', 'TolFun', s.static_tol);
n_eta = floor(s.N / 2);
% initpt = s.alpha * (s.gammaB / (s.gammaA + s.gammaB));
% eta_maxq = fsolve(@(eta) find_argqmax(eta, grid, vA, vB, vAp_c(2:end-1), vBp_c(2:end-1), s), initpt, options);
% eta0 = interp1(grid.eta, grid.eta, eta_maxq, 'nearest');
% n_eta = find(grid.eta == eta0); % location of the initial starting eta, which maximizes q
% if eta0 == eta_maxq
%     grid.eta(n_eta) = eta0 - (eta0 - grid.eta(n_eta-1))/3;
% elseif eta0 > eta_maxq
%     n_eta = n_eta - 1;
% end

% Find minimum capital allocation before specialization stops
s.max_cap = 1 / (1 + s.a / s.a_ * (1 - s.alpha) / s.alpha);

%% Solve from the left first
% Make guess for capital allocation
psiAa0 = s.A_lvg * grid.eta(1);

%% Solve using backward differences
% start with region 1, backward diff
s.side = 0;
for j = 1:n_eta
    eqc = @(psi) eqcond(grid.eta(j), psi, vderivs, grid, s, 1);
    [psiAa1, ~, flag, ~] = fsolve(eqc, psiAa0, options); % options defined above
    [~, grid] = eqcond(grid.eta(j), psiAa1, vderivs, grid, s, 1);
    if grid.psiBa(j) < 0
        break;
    end
    psiAa0 = grid.psiAa(j);
    grid.flag(j) = flag;
end
new_start = j;

% start region 2, backward diff
for j = new_start:n_eta
    eqc = @(psi) eqcond(grid.eta(j), psi, vderivs, grid, s, 2);
    [psiAa1, ~, flag, ~] = fsolve(eqc, psiAa0, options);
    [~, grid] = eqcond(grid.eta(j), psiAa1, vderivs, grid, s, 2);
    grid.flag(j) = flag;
    psiAa0 = grid.psiAa(j);
end

%% Solve using forward differences from right side
s.side = 1;

% make guess
psiAa0 = 1 - (s.B_lvg * (1 - grid.eta(end)));

% start region 3, forward diff
for j = length(grid.eta):-1:(n_eta+1)
    eqc = @(psi) eqcond(grid.eta(j), psi, vderivs, grid, s, 3);
    [psiAa1, ~, flag, ~] = fsolve(eqc, psiAa0, options);
    [~, grid] = eqcond(grid.eta(j), psiAa1, vderivs, grid, s, 3);
    if grid.psiAb(j) < 0
        break;
    end
    psiAa0 = grid.psiAa(j);
    grid.flag(j) = flag;
end
new_start = j;

% start region 2, forward diff
for j = new_start:-1:(n_eta+1)
    eqc = @(psi) eqcond(grid.eta(j), psi, vderivs, grid, s, 2);
    [psiAa1, ~, flag, ~] = fsolve(eqc, psiAa0, options);
    [~, grid] = eqcond(grid.eta(j), psiAa1, vderivs, grid, s, 2);
    psiAa0 = grid.psiAa(j);
    grid.flag(j) = flag;
end

clear vderivs; % Don't need to track these values anymore so clear from memory

%% Fill Grid values
if lastiter == 1
    grid.iota         = ((s.a .* grid.psiAa + s.a_ .* grid.psiBa).^s.alpha .* ...
						(s.a_ .* grid.psiAb + s.a .* grid.psiBb).^(1 - s.alpha) - ...
						s.rho .* grid.q);
    grid.rF           = ((s.a .* grid.Pa - grid.iota) ./ grid.q + grid.muK) - ...
						s.gammaA .* grid.psiA ./ grid.eta .* s.sigA.^2 + (s.gammaA-1) .* ...
                        (s.sigA .* grid.sig_xiA);
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

function [F, grid] = eqcond(eta, psi, vderivs, grid, s, region)
% This function holds the implicit ODE that characterizes
% equilibrium in terms of the capital price q and its first
% derivative qp. This function also computes dynamics
%
% Written by William Chen, Gregory Phelan Mar. 2019

%% Compute capital allocations
psiAa = psi;

% Smoothly interpolate to get vA, vB
vA = interp1(vderivs.eta, vderivs.vA, eta, 'spline', 'extrap');
vB = interp1(vderivs.eta, vderivs.vB, eta, 'spline', 'extrap');
q  = 1/ s.lambda;
xi = vA / eta / q;
% zeta = vB / (1 - eta) / q;

% Find PsiAs, PsiBs, first assuming specialization.  This assumes Cobb-Douglas
% with even weights, otherwise you need to use fsolve.
lvgcons = 0; % indicator for whether an agent is leverage constrained
if region == 2
    % specialization
    if psiAa / eta - 1 > s.LA
        lvgcons = 1;
        psiAa = (1 + s.LA) * eta;
        psiBb = 1 - psiAa;
    elseif (1 - psiAa) / (1 - eta) - 1 > s.LB
        lvgcons = 2;
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
    psiBa = 0; psiAb = 0;
end
if region == 1
    % B doesn't specialize
    prat = s.a / s.a_;
    psiAb = 0;
    psiBa = (1 - (1 + prat) * psiAa) / 2;
    if psiAa / eta - 1 > s.LA
        lvgcons = 1;
        psiAa = (1 + s.LA) * eta;
        psiBa = (1 - (1 + prat) * (1 + s.LA) * eta) / 2;
%     elseif (1 - psiAa) / (1 - eta) - 1 > s.LB
%         lvgcons = 1;
%         psiAa = 1 - (1 + s.LB) * (1 - eta);
%         psiBa = (1 - (1 + prat) * psiAa) / 2;
    end
    psiBb = 1 - psiAa - psiBa;
elseif region == 3
    % A doesn't specialize
    prat = s.a / s.a_;
    psiBa = 0;
    % psiAb = (1 - (1 + prat) * (1 - psiAa - psiAb)) / 2;
    % psiAb = (1 - (1 + prat) * (1 - psiAa)) / 2 + (1+prat) * psiAb / 2;
    psiAb = (1 - (1 + prat) * (1 - psiAa)) / 2 / (1 - (1 + prat) / 2);
    if (1 - psiAa) / (1 - eta) - 1 > s.LB
        lvgcons = 2;
        psiBb = (1 - eta) * (1 + s.LB);
        psiAb = (1 - (1 + prat) * (1 - eta) * (1 + s.LB)) / 2;
        psiAa = 1 - psiBb - psiAb;
    else
        psiBb = 1 - psiAa - psiAb;
    end
end

% Compute output and prices
ya = s.a * psiAa + s.a_ * psiBa;
yb = s.a * psiBb + s.a_ * psiAb;
Y  = ya^s.alpha * yb^(1 - s.alpha);
Pa   = s.alpha * (Y / ya); Pb = (1 - s.alpha) * (Y / yb);
psiA = psiAa + psiAb;
psiB = psiBb + psiBa;

%% Compute volatilities
% Compute capital price volatilities
sig_qA = 0;
sig_qB = 0;
% sig_q = sqrt(sig_qA^2 + sig_qB^2);

% Compute eta volatilities
sig_etaA = psiA * (1 / eta - 1) * s.sigA;
sig_etaB = - psiB * s.sigB;
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
    sig_xiA   = sig_vAA - sig_etaA;
    sig_xiB   = sig_vAB - sig_etaB;
    sig_zetaA = sig_vBA + eta / (1-eta) * sig_etaA;
    sig_zetaB = sig_vBB + eta / (1-eta) * sig_etaB;

    % Compute risk premia
    varsig_Aa    = psiA / eta * s.sigA^2;
    varsig_xiA   = s.sigA * sig_xiA;
    varsig_Bb    = psiB / (1-eta) * s.sigB^2;
    varsig_zetaB = s.sigB * sig_zetaB;

    % Compute drift of eta
    vec_K  = [psiA * s.sigA; psiB * s.sigB];
    vec_Aa = [s.sigA; sig_qB];
    % vec_Bb = [sig_qA; s.sigB + sig_qB];
    iota = (Y - s.rho * q);
    if s.psi == 1
        c_nA = s.rho;
    else
        c_nA = xi^(1 - s.psi);
    end
    if lvgcons == 0 || lvgcons == 2
        mu_eta = (s.a * Pa - iota) / q - c_nA + (psiA / eta - 1) * (s.gammaA * varsig_Aa + ...
            (s.gammaA - 1) * varsig_xiA) + vec_K' * vec_K - ...
            (psiA / eta * vec_Aa)' * vec_K + s.tau * (1 - eta) / eta;
    else
        mu_eta = (psiA / eta - 1) * (s.gammaB * varsig_Bb + (s.gammaB - 1) * varsig_zetaB - ...
            s.a * Pb / q) - iota / q + 1 / eta / q * (psiAa * s.a * Pa + ...
            psiAb * s.a_ * Pb) - c_nA + vec_K' * vec_K - ...
            (psiA / eta * vec_Aa)' * vec_K + s.tau * (1 - eta) / eta;
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
returns_diff   = s.a * (Pa - Pb) / q;
risk_prem_diff = s.gammaA * varsig_Aa + (s.gammaA - 1) * varsig_xiA - ...
                    s.gammaB * varsig_Bb - (s.gammaB - 1) * varsig_zetaB;
if lvgcons == 0
    % then use asset pricing difference to pin down q
    F = returns_diff - risk_prem_diff;
elseif lvgcons == 1
    % then pin down by leverage constraint
    F = (psi + psiAb) - (1 + s.LA) * eta;
elseif lvgcons == 2
    % then pin down by leverage constraint
    F = (1 - psi - psiAb) - (1 + s.LB) * (1 - eta);
end

%% Compute grid values
if nargout > 1
    i = find(grid.eta == eta);
    muK   = iota - s.delta;
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
