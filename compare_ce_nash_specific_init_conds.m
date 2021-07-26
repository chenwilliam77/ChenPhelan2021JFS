% This script computes the consumption equivalent welfare change of the
% competitive equilibrium evaluated at specific initial conditions relative
% to the Nash equilibrium evaluated at the same initial conditions.

print_numbers = 0;

% Initial conditions to compare at
eta0s = linspace(0.01, 0.99, 50)';

% Load CE
load('data/ce.mat');
ce_grid = grid;
ce_s = s;
ce_stat = stat;
ce_welf = welf;

% Load Nash
load('data/nash.mat');
nash_grid = grid;
nash_s = s;
nash_stat = stat;
nash_welf = welf;

% Compute consumption equivalents
ce_expVA = zeros(size(eta0s));
ce_expVB = zeros(size(eta0s));
nash_expVA = zeros(size(eta0s));
nash_expVB = zeros(size(eta0s));
nash2ce_VA = zeros(size(eta0s));
nash2ce_VB = zeros(size(eta0s));
for i = 1:numel(eta0s)
    % Interpolate to get welfare values at each eta
    ce_expVA(i) = interp1(ce_grid.eta, ce_welf.VA, eta0s(i), 'linear');
    ce_expVB(i) = interp1(ce_grid.eta, ce_welf.VB, eta0s(i), 'linear');
    nash_expVA(i) = interp1(nash_grid.eta, nash_welf.VA, eta0s(i), 'linear');
    nash_expVB(i) = interp1(nash_grid.eta, nash_welf.VB, eta0s(i), 'linear');

    % Get consumption equivalent welfare change
    nash2ce_VA(i) = get_consum_equiv(nash_expVA(i), ce_expVA(i), 'A', s);
    nash2ce_VB(i) = get_consum_equiv(nash_expVB(i), ce_expVB(i), 'B', s);

    % Display
    if print_numbers == 1
        disp(['When eta_0 = ', num2str(eta0s(i))]);
        disp('The consumption equivalent equivalent welfare change of ');
        disp('the competitive equilibrium relative to Nash for ');
        disp(['the AE and EME is (', num2str(round(100 * nash2ce_VA(i), 3)), ...
            '%, ', num2str(round(100 * nash2ce_VB(i), 3)), '%)']);
    end
end

ce_sss = ce_grid.eta(find(ce_grid.mu_eta <= 0, 1));
figure(1);
plot(eta0s, 100 * nash2ce_VA); hold on
plot(eta0s, 100 * nash2ce_VB); hold on
plot([ce_sss, ce_sss], [-12, 12], 'color', 'black', 'LineStyle', '--'); hold off
ylabel('pp');
xlabel('$\eta$', 'interpreter', 'latex', 'FontSize', 14);
ylim([-12, 12]);

figure(2);
plot(ce_grid.eta, ce_stat.pdf); hold on
plot(nash_grid.eta, nash_stat.pdf); hold on
max_y = max([max(ce_stat.pdf), max(nash_stat.pdf)]);
plot([ce_sss, ce_sss], [0, max_y * 1.025], 'color', 'black', 'LineStyle', '--'); hold off
ylabel('Stationary Density');
xlabel('$\eta$', 'interpreter', 'latex', 'FontSize', 14);
ylim([0, max_y * 1.025]);
