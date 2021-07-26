% This script analyzes output from hpcc/check_zero.m
%
% Written by William Chen, Jul. 2019

clear; close all;
load('data/hpcc/baseline/check_zero.mat');

expVA = zeros(size(welfs));
expVB = zeros(size(welfs));

for i = 1:numel(welfs)
    expVA(i) = welfs{i}.expVA;
    expVB(i) = welfs{i}.expVB;
end

figure(1);
plot(L0(:,1), expVA(:,1));
xlabel('L_A');
ylabel('E[V_A]');

figure(2);
plot(L0(:,1), expVB(:,2));
xlabel('L_B');
ylabel('E[V_B]');
