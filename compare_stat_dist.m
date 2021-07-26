% compare stationary distributions and examine role of
% systemic risk versus drift of eta.

clear; close all;
addpath parameters;
do_computation = 0; % Set to 0 if you have already saved output in the savefile
savefile = 'data/stat_dist.mat';

if do_computation
    load('data/gamA1p05_gamB2.mat'); % converged solution
    baseline_parameters;
    vAfnct1 = @(x) interp1(grid.eta, welf.vA, x, 'pchip', 'extrap');
    vBfnct1 = @(x) interp1(grid.eta, welf.vB, x, 'pchip', 'extrap');
    s.LA = 1e3; s.LB = 1e3;
    [grid_ce, welf_ce, stat_ce] = get_eqm(vAfnct1, vBfnct1, s, 1, 1, 1);
    s.LA = 1;
    [grid1, welf1, stat1] = get_eqm(vAfnct1, vBfnct1, s, 1, 1, 1);
    save('data/stat_dist.mat');
else
    load(savefile);
end

%% Process data
intce = integrand(grid_ce.eta, grid_ce.mu_eta, grid_ce.sig_etaA, grid_ce.sig_etaB);
intL1 = integrand(grid1.eta, grid1.mu_eta, grid1.sig_etaA, grid1.sig_etaB);
intcemuL1sig = integrand(grid1.eta, grid_ce.mu_eta, grid1.sig_etaA, grid1.sig_etaB);
intL1mucesig = integrand(grid1.eta, grid1.mu_eta, grid_ce.sig_etaA, grid_ce.sig_etaB);

% Compute pdf
[logprece, logdexpce, pdfce] = get_pdf(grid_ce.eta, grid_ce.mu_eta, grid_ce.sig_etaA, grid_ce.sig_etaB);
[logpreL1, logdexpL1, pdfL1] = get_pdf(grid1.eta, grid1.mu_eta, grid1.sig_etaA, grid1.sig_etaB);
[logprecemuL1sig, logdexpcemuL1sig, pdfcemuL1sig] = get_pdf(grid1.eta, grid_ce.mu_eta, grid1.sig_etaA, grid1.sig_etaB);
[logpreL1mucesig, logdexpL1mucesig, pdfL1mucesig] = get_pdf(grid1.eta, grid1.mu_eta, grid_ce.sig_etaA, grid_ce.sig_etaB);

D1ce = (pdfce .* (grid_ce.sig_etaA.^2 + grid_ce.sig_etaB.^2) .* grid_ce.eta.^2);
D1L1 = (pdfL1 .* (grid1.sig_etaA.^2 + grid1.sig_etaB.^2) .* grid1.eta.^2);
D1cemuL1sig = (pdfcemuL1sig .* (grid1.sig_etaA.^2 + grid1.sig_etaB.^2) .* grid1.eta.^2);
D1L1mucesig = (pdfL1mucesig .* (grid_ce.sig_etaA.^2 + grid_ce.sig_etaB.^2) .* grid_ce.eta.^2);

% Compute counterfactual drifts of eta
baseline_parameters;
[L1noToT, ceToT, ceToTq] = analyze_mu(grid1, grid_ce.Pa, grid_ce.Pb, grid_ce.q, 1, s);
[cenoToT, L1ToT, L1ToTq] = analyze_mu(grid_ce, grid1.Pa, grid1.Pb, grid1.q, 1e3, s);
[~, ~, L1mu_cetot] = get_pdf(grid1.eta, ceToT, grid1.sig_etaA, grid1.sig_etaB);
[~, ~, cemu_L1tot] = get_pdf(grid_ce.eta, L1ToT, grid_ce.sig_etaA, grid_ce.sig_etaB);

%% Plot data
figure(1)
plot(grid_ce.eta, intL1 ./ intce); hold on
plot(grid_ce.eta, intcemuL1sig ./ intce); hold on
plot(grid_ce.eta, intL1mucesig ./ intce); hold on
plot(grid_ce.eta, ones(size(grid_ce.eta)), '--', 'color', 'k'); hold off
legend({'L_A = 1 vs. CE', 'CE mu, L_A vol vs. CE', 'L_A mu, CE vol vs. CE'});
xlabel('\eta');

figure(2)
plot(grid_ce.eta, intL1 ./ intce); hold on
plot(grid_ce.eta, intcemuL1sig ./ intce); hold on
plot(grid_ce.eta, intL1mucesig ./ intce); hold on
plot(grid_ce.eta, ones(size(grid_ce.eta)), '--', 'color', 'k'); hold off
legend({'L_A = 1 vs. CE', 'CE mu, L_A vol vs. CE', 'L_A mu, CE vol vs. CE'});
ylim([0;10]);
xlabel('\eta');

figure(3);
plot(grid_ce.eta, logpreL1 ./ logprece); hold on
plot(grid_ce.eta, logprecemuL1sig ./ logprece); hold on
plot(grid_ce.eta, logpreL1mucesig ./ logprece); hold on
plot(grid_ce.eta, ones(size(grid_ce.eta)), '--', 'color', 'k'); hold off
legend({'L_A = 1 vs. CE', 'CE mu, L_A vol vs. CE', 'L_A mu, CE vol vs. CE'});
xlabel('\eta');

figure(4);
plot(grid_ce.eta, logdexpL1 ./ logdexpce); hold on
plot(grid_ce.eta, logdexpcemuL1sig ./ logdexpce); hold on
plot(grid_ce.eta, logdexpL1mucesig ./ logdexpce); hold on
plot(grid_ce.eta, ones(size(grid_ce.eta)), '--', 'color', 'k'); hold off
legend({'L_A = 1 vs. CE', 'CE mu, L_A vol vs. CE', 'L_A mu, CE vol vs. CE'});
xlabel('\eta');

figure(5);
plot(grid_ce.eta, D1L1 ./ D1ce); hold on
plot(grid_ce.eta, D1cemuL1sig ./ D1ce); hold on
plot(grid_ce.eta, D1L1mucesig ./ D1ce); hold on
plot(grid_ce.eta, ones(size(grid_ce.eta)), '--', 'color', 'k'); hold off
legend({'L_A = 1 vs. CE', 'CE mu, L_A vol vs. CE', 'L_A mu, CE vol vs. CE'});
xlabel('\eta');

figure(6);
plot(grid_ce.eta, pdfL1 ./ pdfce); hold on
plot(grid_ce.eta, pdfcemuL1sig ./ pdfce); hold on
plot(grid_ce.eta, pdfL1mucesig ./ pdfce); hold on
plot(grid_ce.eta, ones(size(grid_ce.eta)), '--', 'color', 'k'); hold off
legend({'L_A = 1 vs. CE', 'CE mu, L_A vol vs. CE', 'L_A mu, CE vol vs. CE'});
xlabel('\eta');

figure(7);
plot(grid_ce.eta, D1ce); hold on
plot(grid_ce.eta, D1L1); hold on
plot(grid_ce.eta, D1cemuL1sig); hold on
plot(grid_ce.eta, D1L1mucesig); hold off
legend({'CE', 'L_A = 1', 'CE mu, L_A vol', 'L_A mu, CE vol'});
xlabel('\eta');

figure(8);
plot(grid_ce.eta, pdfce); hold on
plot(grid_ce.eta, pdfL1, '-x'); hold on
plot(grid_ce.eta, pdfcemuL1sig, '--', 'LineWidth', 2); hold on
plot(grid_ce.eta, pdfL1mucesig, '--x', 'LineWidth', 2); hold off
legend({'CE', 'L_A = 1', 'CE mu, L_A vol', 'L_A mu, CE vol'});
xlabel('\eta');
ylabel('Stationary Density')

figure(9);
plot(grid_ce.eta, pdfL1mucesig ./ pdfL1); hold on;
plot(grid_ce.eta, ones(size(grid_ce.eta)), '--', 'color', 'k'); hold off
ylim([-1,1]);
xlabel('\eta');

figure(10);
plot(grid_ce.eta, pdfL1mucesig - pdfL1); hold on
plot(grid_ce.eta, zeros(size(grid_ce.eta)), '--', 'color', 'k'); hold off
xlabel('\eta');
title('PDF of L1 mu CE sig minus PDF of L1 case');
ylabel('Difference');

figure(11);
plot(grid_ce.eta, grid_ce.mu_eta .* grid_ce.eta); hold on
plot(grid1.eta, grid1.mu_eta .* grid1.eta, '-x', 'LineWidth', 2); hold on
% plot(grid1.eta, L1noToT .* grid1.eta); hold on
plot(grid1.eta, ceToT .* grid1.eta, '--'); hold on
% plot(grid1.eta, ceToTq .* grid1.eta, '--'); hold off
xlabel('\eta');
% legend({'L_A = 1 Drift', 'L_A = 1 Drift w/out Price Terms', 'L_A = 1 Drift w/ CE Price Terms'});
% legend({'CE Drift', 'L_A = 1 Drift', 'L_A = 1 Drift w/ CE Price Terms', 'L_A = 1 Drift and q w/ CE Price Terms'});
legend({'CE Drift', 'L_A = 1 Drift', 'L_A = 1 Drift w/ CE Price Terms'});

figure(12);
plot(grid_ce.eta, grid_ce.mu_eta .* grid_ce.eta); hold on
plot(grid1.eta, grid1.mu_eta .* grid1.eta, '-x', 'LineWidth', 2); hold on
% plot(grid_ce.eta, cenoToT .* grid_ce.eta); hold on
plot(grid_ce.eta, L1ToT .* grid_ce.eta, '--'); hold on
% plot(grid_ce.eta, L1ToTq .* grid_ce.eta, '--'); hold off
xlabel('\eta');
% legend({'CE Drift', 'CE Drift w/out Price Terms', 'CE Drift w/ L_A = 1 Price Terms'});
% legend({'CE Drift', 'L_A = 1 Drift', 'CE Drift w/ L_A = 1 Price Terms', 'CE Drift and q w/ L_A = 1 Price Terms'});
legend({'CE Drift', 'L_A = 1 Drift', 'CE Drift w/ L_A = 1 Price Terms'});


figure(13);
plot(grid_ce.eta, pdfce); hold on
plot(grid1.eta, pdfL1, '-x'); hold on
plot(grid_ce.eta, cemu_L1tot, '--', 'LineWidth', 2); hold on
plot(grid1.eta, L1mu_cetot, '--x', 'LineWidth', 2); hold off
xlabel('\eta');
legend({'CE', 'L_A = 1', 'CE Drift w/ L_A = 1 Price Terms', 'L_A = 1 Drift w/ CE Price Terms'});

savefig = figure(12);
saveas(savefig, 'figures/mueta_CEdrift_LA1price', 'epsc');
savefig = figure(11);
saveas(savefig, 'figures/mueta_LA1drift_CEprice', 'epsc');
savefig = figure(13);
saveas(savefig, 'figures/statdist_counterfactual_ce_LA1_LBno', 'epsc');
savefig = figure(8);
saveas(savefig, 'figures/statdist_varyCE_LA1_drift_vol', 'epsc');

function out = integrand(eta, mu_eta, sig_etaA, sig_etaB)
out = mu_eta .* eta ./ ((sig_etaA.^2 + sig_etaB.^2) .* eta.^2);
end
function [logpre, logdexp, pdf] = get_pdf(eta, mu_eta, sig_etaA, sig_etaB)
    integ = @(x) interp1(eta, integrand(eta, mu_eta, sig_etaA, sig_etaB), x, 'linear', 'extrap');
    N = length(eta);
    logpre = zeros(size(eta));
    logdexp = zeros(size(eta));
    logpre(1) = 2 * integral(integ, 0, eta(1));
    logdexp(1) = logpre(1) - log((sig_etaA(1).^2 + sig_etaB(1).^2).* eta(1).^2);

    % Compute pdf everywhere else
    for i = 2:N
        logpre(i) = 2 * integral(integ, 0, eta(i));
        logdexp(i) = logpre(i) - log((sig_etaA(i).^2 + sig_etaB(i).^2) .* eta(i).^2);
    end
    for i = 1:100
        dexp = exp(logdexp / (1 + .01 * (i-1)));
        if sum(dexp == Inf) == 0
            break
        end
    end
    dfnct = @(x) interp1(eta, dexp, x, 'linear', 'extrap');

    % Apply proportionality
    pdf = dexp ./ abs(integral(dfnct, 0, 1));
end

function [noToT, otherToT, otherToTq] = analyze_mu(grid, Pa, Pb, q, LA, s)
    k = find(grid.psiA ./ grid.eta - 1 == LA);
    noToT = grid.mu_eta - s.a * grid.Pa ./ grid.q;
    otherToT = noToT + s.a * Pa ./ grid.q;
    otherToTq = noToT + s.a * Pa ./ q;
    if numel(k) > 0
        noToT(k) = grid.mu_eta(k) + (grid.psiA(k) ./ grid.eta(k) - 1) .* ...
            s.a .* grid.Pb(k) ./ grid.q(k) - 1 ./ grid.eta(k) ./ grid.q(k) .* ...
            (grid.psiAa(k) .* grid.Pa(k) .* s.a +  grid.psiAb(k) .* s.a_ .* grid.Pb(k));
        otherToT(k) = noToT(k) - (grid.psiA(k) ./ grid.eta(k) - 1) .* ...
            s.a .* Pb(k) ./ grid.q(k) + 1 ./ grid.eta(k) ./ grid.q(k) .* ...
            (grid.psiAa(k) .* Pa(k) .* s.a +  grid.psiAb(k) .* s.a_ .* Pb(k));
        otherToTq(k) = noToT(k) - (grid.psiA(k) ./ grid.eta(k) - 1) .* ...
            s.a .* Pb(k) ./ q(k) + 1 ./ grid.eta(k) ./ q(k) .* ...
            (grid.psiAa(k) .* Pa(k) .* s.a +  grid.psiAb(k) .* s.a_ .* Pb(k));
    end
end
