function [grid, welf, stat] = get_eqm(vAfnct, vBfnct, s, yes_stat, yes_display, yes_time)

if nargin == 3
   yes_stat = 0;
end
if nargin <= 4
    yes_display = 0;
end
if nargin <= 5
    yes_time = 0;
end
if not(isfield(s, 'Viter'))
    % Default convergence iterations
    s.Viter = 40;
end

% Initialize grid and guess
grid = init_grid(s);
vA = vAfnct(grid.eta);
vB = vBfnct(grid.eta);

welf.vAerr_all = zeros(length(grid.eta), s.Viter);
welf.vBerr_all = zeros(length(grid.eta), s.Viter);
L2_err = zeros(s.Viter, 2);
Linf_err = zeros(s.Viter, 2);
if yes_display == 1
    disp('Solving for equilibrium . . .');
end
if yes_time == 1
    eqm_tim = tic();
end
if isfield(s, 'lvg_approx') && s.lvg_approx == 1
    % We approximate the competitive equilibrium by raising a binding
    % leverage constraint until it stops binding except at the endpoint
    s.LB = s.LB0;
    for il = 1:s.lvg_approx_iter
        if yes_display == 1
            disp(['Iteration ' num2str(il) ' of ' num2str(s.lvg_approx_iter) ', constraint guess = ' num2str(s.LB)]);
        end
        for t = 1:s.Viter
            if mod(t, 10) == 0 && yes_display == 1
                tmpNum = t / s.Viter * 100;
                disp(['Percentage complete: ' num2str(tmpNum) '%']);
            end

            % Compute static step
            if ~isreal(vA) || ~isreal(vB)
                keyboard;
            end
            [grid, s] = staticstep_fd(grid, vA, vB, s, 0);

            % Compute PDE step
            [grid, vA1, vB1] = pdestep(grid, vA, vB, s);
            vA1 = real(vA1); vB1 = real(vB1);

            % Compute errors
            vAerr = reshape((vA1 - vA) ./ vA1, numel(vA1), 1);
            vBerr = reshape((vB1 - vB) ./ vB1, numel(vB1), 1);
            welf.vAerr_all(:,t) = vAerr;
            welf.vBerr_all(:,t) = vBerr;
            L2_err(t,1)   = norm(vAerr); % vA1 - vA ~ vA_t, divide by vA1 to get implied error
            Linf_err(t,1) = max(abs(vAerr));
            L2_err(t,2)   = norm(vBerr);
            Linf_err(t,2) = max(abs(vBerr));


            if max(max(Linf_err(t,:))) < s.V_Linftol
                if yes_display == 1
                    disp('Convergence achieved for value functions');
                end
                [grid, s] = staticstep_fd(grid, vA1, vB1, s, 1);
                vA = vA1;
                vB = vB1;
                break;
            end

            % Update value functions
            vA = vA1;
            vB = vB1;
        end

        lvgB = grid.psiB ./ (1 - grid.eta) - 1;
        if length(find(lvgB == s.LB)) <= 1 - s.lvg_approx_full
            if max(max(Linf_err(t,:))) < s.V_Linftol
                if yes_display == 1
                    disp('Competitive equilibrium successfully computed');
                end
            end
            break;
        end
        if grid.eta(find(lvgB == s.LB, 1)) >= s.end - s.lvg_approx_switchLinc
            s.LB = s.LB + s.lvg_approx_Linc2;
        else
            s.LB = s.LB + s.lvg_approx_Linc1;
        end
    end
    if il == s.lvg_approx_iter && yes_display == 1
        disp('Failed to compute competitive equilibrium');
    end
else
    for t = 1:s.Viter
        if mod(t, 10) == 0 && yes_display == 1
            tmpNum = t / s.Viter * 100;
            disp(['Percentage complete: ' num2str(tmpNum) '%']);
            figure(10000);
            plot(grid.eta, vA); hold on
            figure(10001);
            plot(grid.eta, vB); hold on
        end

        % Compute static step
        if ~isreal(vA) || ~isreal(vB)
            keyboard;
        end
        
        if t < s.Viter
            [grid, s] = staticstep_fd(grid, vA, vB, s, 0);
        else
            % Then we are at the last iteration, and we want to compute
            % equilibrium data
            [grid, s] = staticstep_fd(grid, vA, vB, s, 1);
        end

        % Compute PDE step
        [grid, vA1, vB1] = pdestep(grid, vA, vB, s);
        vA1 = real(vA1); vB1 = real(vB1);

        % Compute errors
        vAerr = reshape((vA1 - vA) ./ vA1, numel(vA1), 1);
        vBerr = reshape((vB1 - vB) ./ vB1, numel(vB1), 1);
        welf.vAerr_all(:,t) = vAerr;
        welf.vBerr_all(:,t) = vBerr;
        L2_err(t,1)   = norm(vAerr); % vA1 - vA ~ vA_t, divide by vA1 to get implied error
        Linf_err(t,1) = max(abs(vAerr));
        L2_err(t,2)   = norm(vBerr);
        Linf_err(t,2) = max(abs(vBerr));
        welf.vAerr = vAerr;
        welf.vBerr = vBerr;
    %     disp(t);
        if max(max(Linf_err(t,:))) < s.V_Linftol
            if yes_display == 1
                disp('Convergence achieved for value functions');
            end
            [grid, s] = staticstep_fd(grid, vA, vB, s, 1);
            vA = vA1;
            vB = vB1;
            welf.V_converge = 1;
            break;
        end
        if t == s.Viter
            if yes_display == 1
                disp('Convergence not achieved for value functions');
            end
            welf.V_converge = 0;
        end

        % Update value functions
        vA = vA1;
        vB = vB1;
    end
end
welf.vA = vA;
welf.vB = vB;
welf.VA = vA.^(1-s.gammaA) ./ (1-s.gammaA);
welf.VB = vB.^(1-s.gammaB) ./ (1-s.gammaB);
welf.L2_err   = L2_err(1:t,:);
welf.Linf_err = Linf_err(1:t,:);
grid.xi   = welf.vA ./ grid.eta ./ grid.q;
grid.zeta = welf.vB ./ (1 - grid.eta) ./ grid.q;
if yes_time == 1
    toc(eqm_tim);
end

if yes_stat == 1
    if yes_display == 1
       disp('Computing stationary density . . .');
    end
    if yes_time == 1
       stat_tim = tic();
    end
    [stat.pdf, stat.cdf] = ergocalc(grid);
    if yes_display == 1
       disp('Computing ex-ante welfare . . .');
    end
    integ      = @(x) interp1(grid.eta, stat.pdf .* welf.VA, x, 'linear', 'extrap');
    welf.expVA = integral(integ, grid.eta(1), grid.eta(end));
    integ      = @(x) interp1(grid.eta, stat.pdf .* welf.VB, x, 'linear', 'extrap');
    welf.expVB = integral(integ, grid.eta(1), grid.eta(end));
    if yes_time == 1
       toc(stat_tim);
    end
end

% Compute calibration statistics
grid.avg_rF = integral(@(x) interp1(grid.eta, grid.rF, x, 'pchip', 'extrap'), ...
        grid.eta(1), grid.eta(end));
stoch_ss = find(grid.mu_eta <= 0, 1);
grid.ca_gdp = -((grid.psiA - grid.eta) .* grid.q) ./ ...
              (s.a * grid.Pa .* grid.psiAa + s.a_ * grid.Pb .* grid.psiAb);
grid.ca_gdp_sss = grid.ca_gdp(stoch_ss);
grid.invst_gdp = (s.lambda * grid.q - 1) ./ s.kappa ./ ...
                ((s.a .* grid.psiAa + s.a_ * grid.psiBa).^s.alpha .* ...
                (s.a_ .* grid.psiAb + s.a * grid.psiBb).^(1 - s.alpha));
grid.avg_invst_gdp = integral(@(x) interp1(grid.eta, grid.invst_gdp, x, 'pchip', 'extrap'), ...
                    grid.eta(1), grid.eta(end));
grid.output = ((s.a .* grid.psiAa + s.a_ * grid.psiBa).^s.alpha .* ...
                (s.a_ .* grid.psiAb + s.a * grid.psiBb).^(1 - s.alpha));
grid.consum = grid.output - grid.iota;
grid.iota = grid.output - s.rho;
grid.growth = grid.iota - s.delta;
grid.avg_growth = integral(@(x) interp1(grid.eta, grid.growth, x, 'pchip', 'extrap'), ...
                    grid.eta(1), grid.eta(end));
varsigAa  = grid.psiA ./ grid.eta .* ((s.sigA + grid.sig_qA).^2 + grid.sig_qB.^2);
varsigxiA = grid.sig_xiA .* (s.sigA + grid.sig_qA) + grid.sig_xiB .* grid.sig_qB;
varsigBb  = grid.psiB ./ (1 - grid.eta) .* (grid.sig_qA.^2 + (s.sigB + grid.sig_qB).^2);
varsigzetaB = grid.sig_zetaA .* grid.sig_qA + grid.sig_zetaB .* (s.sigB + grid.sig_qB);
grid.ex_retA = s.gammaA .* varsigAa + (s.gammaA - 1) .* varsigxiA;
grid.ex_retB = s.gammaB .* varsigBb + (s.gammaB - 1) .* varsigzetaB;
grid.avg_ex_retA = integral(@(x) interp1(grid.eta, grid.ex_retA, x, 'pchip', 'extrap'), ...
                    grid.eta(1), grid.eta(end));
grid.avg_ex_retB = integral(@(x) interp1(grid.eta, grid.ex_retB, x, 'pchip', 'extrap'), ...
                    grid.eta(1), grid.eta(end));
grid.avg_diff_ex_ret = integral(@(x) interp1(grid.eta, abs(grid.ex_retA - grid.ex_retB), x, 'pchip', 'extrap'), ...
                    grid.eta(1), grid.eta(end));
grid.lvgA = grid.psiA ./ grid.eta - 1;
grid.lvgB = grid.psiB ./ (1 - grid.eta) - 1;
last_i = find(grid.lvgA < 0, 1);
grid.avg_lvgA = integral(@(x) interp1(grid.eta, grid.lvgA, x, 'pchip', 'extrap'), ...
                    grid.eta(1), grid.eta(last_i-1)) / ...
                    (grid.eta(last_i-1) - grid.eta(1));
grid.avg_lvgB = integral(@(x) interp1(grid.eta, grid.lvgB, x, 'pchip', 'extrap'), ...
                    grid.eta(last_i), grid.eta(end)) / ...
                    (grid.eta(end) - grid.eta(last_i));
grid.sharpeA = grid.ex_retA ./ sqrt(varsigAa);
grid.sharpeB = grid.ex_retB ./ sqrt(varsigBb);
grid.avg_sharpeA = integral(@(x) interp1(grid.eta, grid.sharpeA, x, 'pchip', 'extrap'), ...
                    grid.eta(1), grid.eta(end));
grid.avg_sharpeB = integral(@(x) interp1(grid.eta, grid.sharpeB, x, 'pchip', 'extrap'), ...
                    grid.eta(1), grid.eta(end));

if yes_display == 1
    disp(['Average risk-free interest rate: ' num2str(grid.avg_rF)]);
    disp(['Current account deficit at stochastic steady state: ' num2str(grid.ca_gdp_sss)]);
    disp(['Average growth rate: ' num2str(grid.avg_growth)]);
    disp(['Average investment to output ratio: ' num2str(grid.avg_invst_gdp)]);
    disp(['Average excess returns for A: ' num2str(grid.avg_ex_retA)]);
    disp(['Average excess returns for B: ' num2str(grid.avg_ex_retB)]);
    disp(['Average difference in excess returns: ' num2str(grid.avg_diff_ex_ret)]);
    disp(['Average leverage for A: ' num2str(grid.avg_lvgA)]);
    disp(['Average leverage for B: ' num2str(grid.avg_lvgB)]);
    disp(['Average Sharpe ratio for A: ' num2str(grid.avg_sharpeA)]);
    disp(['Average Sharpe ratio for B: ' num2str(grid.avg_sharpeB)]);
    disp(['Expected welfare for A: ' num2str(welf.expVA)]);
    disp(['Expected welfare for B: ' num2str(welf.expVB)]);
end
end

