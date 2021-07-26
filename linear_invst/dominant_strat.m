% This script illustrates the effects of removing endogenous asset price volatility
% on the expected welfare for country A.
%
% Written by William Chen, Jun. 2019

clear; close all;
load('data/gamA1p05_gamB2'); % converged solution
linear_invst_baseline_parameters;
vAfnct1 = @(x) interp1(grid.eta, welf.vA, x, 'pchip', 'extrap');
vBfnct1 = @(x) interp1(grid.eta, welf.vB, x, 'pchip', 'extrap');

L0 = [1e3; linspace(1, 0, 11)'];
% expVA = zeros(size(L0));
% s.LB = 1e3;
% ss = cell(size(L0));
% for i = 1:length(L0)
%     ss{i} = s;
%     ss{i}.LA = L0(i);
% end
% for i = 1:length(L0)
%     [~, welf, ~] = get_eqm(vAfnct1, vBfnct1, ss{i}, 1, 0, 0);
%     expVA(i) = welf.expVA;
% end
%
% figure(1);
% plot(L0(2:end), expVA(2:end)); hold on
% plot(L0(2:end), ones(size(L0(2:end))) * expVA(1), '--'); hold on
% xlabel('L_A');
% title('Expected Welfare for Country A');
% legend({'Binding L_A', 'CE'}, 'location', 'northeast', 'FontSize', 16);

expVA = zeros(size(L0));
s.LB = 0;
ss = cell(size(L0));
for i = 1:length(L0)
    ss{i} = s;
    ss{i}.LA = L0(i);
end
for i = 1:length(L0)
    [~, welf, ~] = get_eqm(vAfnct1, vBfnct1, ss{i}, 1, 0, 0);
    expVA(i) = welf.expVA;
end

figure(2);
plot(L0(2:end), expVA(2:end)); hold on
plot(L0(2:end), ones(size(L0(2:end))) * expVA(1), '--'); hold on
xlabel('L_A');
% title('Expected Welfare for Country A');
legend({'Binding L_A', 'CE'}, 'location', 'northeast', 'FontSize', 16);
savefig = figure(2);
saveas(savefig, 'figures/expVA_varyLA_LB0', 'epsc');

% expVB = zeros(size(L0));
% s.LA = 0;
% ss = cell(size(L0));
% for i = 1:length(L0)
%     ss{i} = s;
%     ss{i}.LB = L0(i);
% end
% for i = 1:length(L0)
%     [~, welf, ~] = get_eqm(vAfnct1, vBfnct1, ss{i}, 1, 0, 0);
%     expVB(i) = welf.expVB;
% end
%
% figure(3);
% plot(L0(2:end), expVB(2:end)); hold on
% plot(L0(2:end), ones(size(L0(2:end))) * expVB(1), '--'); hold on
% xlabel('L_B');
% title('Expected Welfare for Country B');
% legend({'Binding L_B', 'CE'}, 'location', 'northeast', 'FontSize', 16);
%
%
