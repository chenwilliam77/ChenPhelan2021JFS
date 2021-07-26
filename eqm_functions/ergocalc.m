function [pdf, cdf] = ergocalc(grid)
% This function calculates the ergodic (stationary) distribution
% using proportionality and numerical integration. The equation
% is derived by using Kolmogorov's forward equation, boundary conditions
% for a non-degenerate stationary distribution, and the integrating factor method.
%
% Written by William Chen, Mar. 2019

% Set up integrand
N = grid.dim1;
integrand = @(x) interp1(grid.eta, (grid.mu_eta .* grid.eta) ./ ...
	((grid.sig_etaA.^2 + grid.sig_etaB.^2) .* ...
	grid.eta.^2), x, 'linear', 'extrap');

% Compute pdf
logdexp = zeros(N,1);

% Compute pdf at eta(1)
logdexp(1) = 2 * integral(integrand, 0, grid.eta(1)) - log((grid.sig_etaA(1).^2 + grid.sig_etaB(1).^2).* grid.eta(1).^2);

% Compute pdf everywhere else
for i = 2:N
    logdexp(i) = 2 * integral(integrand, 0, grid.eta(i)) - log((grid.sig_etaA(i).^2 + grid.sig_etaB(i).^2).* grid.eta(i).^2);
end
for i = 1:100000
    dexp = exp(logdexp / (1 + .01 * (i-1)));
    if sum(dexp == Inf) == 0
        try
			normalize_C = integral(@(x) interp1(grid.eta, dexp, x, 'linear', 'extrap'), 0, 1);
            if normalize_C ~= Inf && ~isnan(normalize_C)
		        break
            end
        catch
            continue;
        end
    end
end
% dexp = zeros(N,1);
% dexp(1) = 1 ./ ((grid.sig_etaA(1).^2 + grid.sig_etaB(1).^2).* grid.eta(1).^2) ...
%         .* exp(2 * integral(integrand, 0, grid.eta(1)));
% for i = 2:N
%     dexp(i) = 1 ./ ((grid.sig_etaA(i).^2 + grid.sig_etaB(i).^2).* grid.eta(i).^2) ...
%         .* exp(2 * integral(integrand, 0, grid.eta(i)));
% end
dfnct = @(x) interp1(grid.eta, dexp, x, 'linear', 'extrap');

% Apply proportionality
pdf = dexp ./ abs(integral(dfnct, 0, 1));

% Compute CDF
cdf = zeros(N, 1);
integ = @(x) interp1(grid.eta, pdf, x, 'linear', 'extrap');
for i = 1:N
    cdf(i) = integral(integ, 0, grid.eta(i));
end

end
