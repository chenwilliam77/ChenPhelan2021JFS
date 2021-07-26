function out = square_root_interval_to_real_line(x, a, b, c)
% out = square_root_interval_to_real_line(x, a, b, c)
%
% This function transforms x from the interval [a, b]
% to the real line using a square root transformation
% See line 900 of https://github.com/FRBNY-DSGE/ModelConstructors.jl/blob/main/src/parameters.jl
%
% Written by William Chen, May 2021
    if nargin < 4
        c = 1;
    end

    cx  = 2 * (x - (a + b) / 2) / (b - a);
    out = (1 / c) * cx / sqrt(1 - cx^2);
end
