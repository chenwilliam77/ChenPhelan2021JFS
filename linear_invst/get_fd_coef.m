function out = get_fd_coef(deriv, acc, gridpts, evalpts)
out = cell2mat(arrayfun(@(x0) wrap_fd(deriv, acc, x0, gridpts), evalpts, 'UniformOutput', false));
end

function out = weights(deriv, x0, gridpt, all)
% This function implements Fornberg's algorithm
% for computing finite difference coefficients
% of arbitrary order and accuracy on non-uniform grids.
% deriv: order of derivative
% acc: order of accuracy
% x0: point at which the finite difference is evaluated, not required to be an element of gridpt
% gridpt: coordinates of non-uniform grid points, e.g. Chebyshev nodes, [1, 1.2, 1.5],
% all: 1 if all coefficients, including lower derivatives than
%      specified are desired.
%
% Written by William Chen and Michael Curran, Apr. 2019

    if nargin < 4
        all = 0;
    end
    if deriv < 0
		error('Order of derivative canot be negative');
    end
    n = length(gridpt); c = zeros(deriv + 1, n); % c is coefficients
    c1 = 1; c4 = gridpt(1) - x0; c(1,1)=1;
    for i = 2:n
       mn = min(i, deriv + 1); c2 = 1; c5 = c4; c4 = gridpt(i) - x0;
       for j = 1:i-1
          c3 = gridpt(i) - gridpt(j);  c2 = c2 * c3;
          if j == i-1 
             c(2:mn,i) = c1*((1:mn-1)' .* c(1:mn-1,i-1) - c5 * c(2:mn,i-1)) / c2;
             c(1,i) = -c1 * c5 *c(1,i-1) / c2;
          end
          c(2:mn,j) =(c4 * c(2:mn,j) - (1:mn-1)' .*c (1:mn-1,j)) / c3;
          c(1,j) = c4 * c(1,j) / c3;
       end
       c1 = c2;
    end

    % Decide what to output
	if all == 1
		out = c;
    else
        out = c(end, :); % assume want requested derivative
	end
end

function out = wrap_fd(deriv, acc, x0, gridpt, all)
    if nargin < 6
        all = 0;
    end
    i = find(gridpt == x0);
    % central difference
    if deriv == 1
        grid_samp = gridpt(i-1:i+acc/2);
    elseif i == 2 && deriv == 2
        grid_samp = gridpt(i-1:i+acc);
    else
        grid_samp = gridpt(i-acc:i+1);
    end
    out = weights(deriv, x0, grid_samp, all);
end