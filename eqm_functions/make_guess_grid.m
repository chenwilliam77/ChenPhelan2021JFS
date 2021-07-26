function [guessA, guessB] = make_guess_grid(eta0, vA, vB)
% This function creates a cell of value functions to be used for
% finding a decent close guess. The function assumes the user
% knows the look up structure to inputs vA and vB.
%
% Written by William Chen, Jun. 2019
    guessA = cell(size(vA));
    guessB = cell(size(vB));
    for i = 1:numel(vA)
        % create grid of guesses
        guessA{i} = @(x) interp1(eta0, vA{i}, x, 'pchip', 'extrap');
        guessB{i} = @(x) interp1(eta0, vB{i}, x, 'pchip', 'extrap');
    end   
end