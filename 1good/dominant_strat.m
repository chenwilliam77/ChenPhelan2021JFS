% This script illustrates the effects of the terms of trade hedge
% on the expected welfare for country A.
%
% Written by William Chen, Jun. 2019

clear; close all;
load('data/gamA1p05_gamB2'); % converged solution
one_good_baseline_parameters;
vAfnct1 = @(x) interp1(grid.eta, welf.vA, x, 'pchip', 'extrap');
vBfnct1 = @(x) interp1(grid.eta, welf.vB, x, 'pchip', 'extrap');

L0 = [1e3; linspace(1, 0, 10)'];
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
% legend({'Binding L_A', 'CE'}, 'location', 'southeast', 'FontSize', 16);

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
legend({'Binding L_A', 'CE'}, 'location', 'east', 'FontSize', 16);
savefig = figure(2);
saveas(savefig, 'figures/expVA_varyLA_LB0.eps', 'epsc');
