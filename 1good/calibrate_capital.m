% This script chooses a, lambda, and kappa to calibrate the model
%
% Written by William Chen, May 2019

one_good_baseline_parameters;
targ1 = .025;
targ2 = .17;
answer_max = fsolve(@(x) calibrate_max(x, [targ1; targ2], s), [1;1]); % .16 based on US data
% answer_min = fsolve(@(x) calibrate_min(x, [targ1; targ2], s), [1;1]);
% answer_a   = fsolve(@(x) calibrate_a(x, [targ1; targ2]], s), [1;1]);
% answer_all = fsolve(@(x) calibrate_all(x, [targ1; targ2], s), [1;1;1]);
disp(answer_max');
% disp(answer_min');
% disp(answer_a');
% disp(answer_all');

function out = calibrate_max(x, targ, s)
    lambda = x(1);
    kappa  = x(2);
    qmax = (kappa * s.a + 1) / (s.rho * kappa + lambda);
    out = zeros(2,1);
    out(1) = lambda / kappa * log(lambda * qmax) - s.delta - targ(1);
    out(2) = (lambda * qmax - 1) / (kappa * s.a) - targ(2);
end
function out = calibrate_min(x, targ, s)
    lambda = x(1);
    kappa  = x(2);
    qmin = (kappa * sqrt(s.a * s.a_) / 2 + 1) / (s.rho * kappa + lambda);
    out = zeros(2,1);
    out(1) = lambda / kappa * log(lambda * qmin) - s.delta - targ(1);
    out(2) = (lambda * qmin - 1) / (kappa * s.a / 2) - targ(2);
end

function out = calibrate_a(x, targ, s)
    a = x(1);
    kappa  = x(2);
    qmax = (kappa * a / 2 + 1) / (s.rho * kappa + s.lambda);
    out = zeros(2,1);
    out(1) = s.lambda / kappa * log(s.lambda * qmax) - s.delta - targ(1);
    out(2) = (s.lambda * qmax - 1) / (kappa * a / 2) - targ(2);
end

function out = calibrate_all(x, targ, s)
    lambda = x(1);
    kappa  = x(2);
    a = x(3);
    qmax = (kappa * a / 2 + 1) / (s.rho * kappa + lambda);
    out = zeros(2,1);
    out(1) = lambda / kappa * log(lambda * qmax) - s.delta - targ(1);
    out(2) = (lambda * qmax - 1) / (kappa * a / 2) - targ(2);
end
