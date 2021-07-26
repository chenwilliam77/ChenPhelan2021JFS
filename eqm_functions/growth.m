function out = growth(q, s)
% This function returns the capital growth rate depending on q
%
% Written by William Chen, Gregory Phelan 2017

% Phi(iota) = log(kappa * iota + 1), w/ Phi'(iota) = 1 / q
iota = 1 / s.kappa .* (s.lambda .* q - 1);
out = s.lambda / s.kappa .* log(s.kappa .* iota + 1);

% qmax = max(qvec);
% qmin = min(qvec);
% frac = par.frac; %this tells us where to make the cutoff
% qmin=min(qvec);
% qmax=max(qvec);
% qbar=frac*qmin+(1-frac)*qmax;
% zeta = par.zeta;
% curvature = par.curv;
% 
% if q<qbar,
%     gq=0;
% elseif curvature == 0,
%     gq = par.g;
% %     g = par.g;
% %     g_ = par.g_;
% %     gq = abs((g-g_)*(q-qmin)/(qmax-qmin));
% 
% else
%    g = par.g;
%    gq = g*((q-qmin)/(qmax-qmin))^zeta;
% end
end
