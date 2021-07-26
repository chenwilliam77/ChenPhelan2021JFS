% compare stationary distributions and examine role of
% systemic risk versus drift of eta.

clear; close all;
load('stat_dist.mat');

%% Process data
intce = integrand(grid_ce.eta, grid_ce.mu_eta, grid_ce.sig_etaA, grid_ce.sig_etaB);
intL1 = integrand(grid_1.eta, grid_1.mu_eta, grid_1.sig_etaA, grid_1.sig_etaB);
intcemuL1sig = integrand(grid_1.eta, grid_ce.mu_eta, grid_1.sig_etaA, grid_1.sig_etaB);
intL1mucesig = integrand(grid_1.eta, grid_1.mu_eta, grid_ce.sig_etaA, grid_ce.sig_etaB);
intLp3 = integrand(grid_p3.eta, grid_p3.mu_eta, grid_p3.sig_etaA, grid_p3.sig_etaB);
intcemuLp3sig = integrand(grid_p3.eta, grid_ce.mu_eta, grid_p3.sig_etaA, grid_p3.sig_etaB);
intLp3mucesig = integrand(grid_p3.eta, grid_p3.mu_eta, grid_ce.sig_etaA, grid_ce.sig_etaB);

% Compute pdf
[logprece, logdexpce, pdfce] = get_pdf(grid_ce.eta, grid_ce.mu_eta, grid_ce.sig_etaA, grid_ce.sig_etaB);
[logpreL1, logdexpL1, pdfL1] = get_pdf(grid_1.eta, grid_1.mu_eta, grid_1.sig_etaA, grid_1.sig_etaB);
[logpreLp1, logdexpLp1, pdfLp1] = get_pdf(grid_p1.eta, grid_p1.mu_eta, grid_p1.sig_etaA, grid_p1.sig_etaB);
[logprecemuL1sig, logdexpcemuL1sig, pdfcemuL1sig] = get_pdf(grid_1.eta, grid_ce.mu_eta, grid_1.sig_etaA, grid_1.sig_etaB);
[logpreL1mucesig, logdexpL1mucesig, pdfL1mucesig] = get_pdf(grid_1.eta, grid_1.mu_eta, grid_ce.sig_etaA, grid_ce.sig_etaB);
[logpreLp3, logdexpLp3, pdfLp3] = get_pdf(grid_p3.eta, grid_p3.mu_eta, grid_p3.sig_etaA, grid_p3.sig_etaB);
[logprecemuLp3sig, logdexpcemuLp3sig, pdfcemuLp3sig] = get_pdf(grid_p3.eta, grid_ce.mu_eta, grid_p3.sig_etaA, grid_p3.sig_etaB);
[logpreLp3mucesig, logdexpLp3mucesig, pdfLp3mucesig] = get_pdf(grid_p3.eta, grid_p3.mu_eta, grid_ce.sig_etaA, grid_ce.sig_etaB);

D1ce = (pdfce .* (grid_ce.sig_etaA.^2 + grid_ce.sig_etaB.^2) .* grid_ce.eta.^2);
D1L1 = (pdfL1 .* (grid_1.sig_etaA.^2 + grid_1.sig_etaB.^2) .* grid_1.eta.^2);
D1cemuL1sig = (pdfcemuL1sig .* (grid_1.sig_etaA.^2 + grid_1.sig_etaB.^2) .* grid_1.eta.^2);
D1L1mucesig = (pdfL1mucesig .* (grid_ce.sig_etaA.^2 + grid_ce.sig_etaB.^2) .* grid_ce.eta.^2);

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
plot(grid_ce.eta, pdfL1); hold on
plot(grid_ce.eta, pdfcemuL1sig, '--'); hold on
plot(grid_ce.eta, pdfL1mucesig, '--'); hold off
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
plot(grid_ce.eta, pdfce); hold on
plot(grid_ce.eta, pdfLp3); hold on
plot(grid_ce.eta, pdfcemuLp3sig, '--'); hold on
plot(grid_ce.eta, pdfLp3mucesig, '--'); hold off
legend({'CE', 'L_A = 0.3', 'CE mu, L_A vol', 'L_A mu, CE vol'});
xlabel('\eta');
ylabel('Stationary Density')

figure(12);
plot(grid_ce.eta, pdfce); hold on
plot(grid_1.eta, pdfL1); hold on
plot(grid_1.eta, pdfLp3); hold on
plot(grid_1.eta, pdfLp1); hold on
legend({'CE', 'L_A = 1', 'L_A = 0.3', 'L_A = 0.1'});
xlabel('\eta');
title('Stationary Density')

figure(13);
plot(grid_ce.eta, grid_ce.mu_eta .* grid_ce.eta); hold on
plot(grid_1.eta, grid_1.mu_eta .* grid_1.eta); hold on
plot(grid_p3.eta, grid_p3.mu_eta .* grid_p3.eta); hold on
plot(grid_p1.eta, grid_p1.mu_eta .* grid_p1.eta); hold on
legend({'CE', 'L_A = 1', 'L_A = 0.3', 'L_A = 0.1'});
xlabel('\eta');
title('\mu_\eta\eta')

figure(14);
plot(grid_ce.eta, sqrt(grid_ce.sig_etaA.^2 + grid_ce.sig_etaB.^2) .* grid_ce.eta); hold on
plot(grid_1.eta, sqrt(grid_1.sig_etaA.^2 + grid_1.sig_etaB.^2) .* grid_1.eta); hold on
plot(grid_p3.eta, sqrt(grid_p3.sig_etaA.^2 + grid_p3.sig_etaB.^2) .* grid_p3.eta); hold on
plot(grid_p1.eta, sqrt(grid_p1.sig_etaA.^2 + grid_p1.sig_etaB.^2) .* grid_p1.eta); hold on
legend({'CE', 'L_A = 1', 'L_A = 0.3', 'L_A = 0.1'});
xlabel('\eta');
title('Volatility of \eta')


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
    for i = 1:1000
        dexp = exp(logdexp / (1 + .01 * (i-1)));
        if sum(dexp == Inf) == 0
            break
        end
    end
    dfnct = @(x) interp1(eta, dexp, x, 'linear', 'extrap');

    % Apply proportionality
    pdf = dexp ./ abs(integral(dfnct, 0, 1));
end
