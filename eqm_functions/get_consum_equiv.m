function phi = get_consum_equiv(v0, v1, country, s)
% This function computes the percentage change in welfare
% in consumption equivalent units, where v0 is the value function
% in the original equilibrium and v1 is the value function
% in the alternative equilibrium. The output gives both gains
% and costs, relative to the initial equilibrium.
%
% Written by William Chen, Jun. 2019
    if country == 'A'
        if v1 > v0
            phi = (v1 ./ v0).^(1 / (1 - s.gammaA)) - 1;
        else
            phi = -(1 - (v1 ./ v0).^(1 / (1 - s.gammaA)));
        end
    elseif country == 'B'
        if v1 > v0
            phi = (v1 / v0).^(1 / (1 - s.gammaB)) - 1;
        else
            phi = -(1 - (v1 / v0).^(1 / (1 - s.gammaB)));
        end
    end
end
