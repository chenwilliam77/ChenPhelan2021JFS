function out = square_root_real_line_to_interval(x, a, b, c)
% out = square_root_real_line_to_interval(x, a, b, c)
%
% This function transforms x from the real line to the
% interval [a, b] using a square root transformation
% See line 765 of https://github.com/FRBNY-DSGE/ModelConstructors.jl/blob/main/src/parameters.jl
%
% Written by William Chen, May 2021
    if nargin < 4
        c = 1;
    end

    out = (a + b) / 2 + ((b - a) / 2) * (c * x) / sqrt(1 + c^2 * x^2);
end
