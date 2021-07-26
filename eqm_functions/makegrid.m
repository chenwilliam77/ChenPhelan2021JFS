function out = makegrid(x0, L, N)
if N <= 1
    error('N needs to be greater than 1');
end
out = x0 + 1 / 2 * L * (1 - cos(pi * (0:(N-1))' / (N-1)));
end