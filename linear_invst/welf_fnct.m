function fp = welf_fnct(eta,f,etaout,q,mu_eta,sigma_etaa,sigma_etab,psiAa,psiAb,par)

% This function solves the second order ODE for HA.
% f = [HA, y], y = HA', fp = [HA',y']
%
% Written by William Chen, Gregory Phelan 2017

rho = par.rho;
sigmaA = par.sigmaA;
sigmaB = par.sigmaB;

interp = 'linear';
mu_eta1 = interp1(etaout,mu_eta,eta,interp);
q1 = interp1(etaout,q,eta,interp);
sigma_etaa1 = interp1(etaout,sigma_etaa,eta,interp);
sigma_etab1 = interp1(etaout,sigma_etab,eta,interp);
psiAa1 = interp1(etaout,psiAa,eta,interp);
psiAb1 = interp1(etaout,psiAb,eta,interp);
fp = zeros(2,1);
% Phi = (q1-1)/kappa;

%muK = growth_q(q1,q,par);
%muK = growth_q1(q1,q,par);
muK = par.g;
%muK = par.g;
psiA1 = psiAa1 + psiAb1;
psiB1 = 1-psiA1;
fp(1) = f(2);
fp(2) = 2/eta^2/(sigma_etaa1^2+sigma_etab1^2)*(rho*f(1)-...
    log(rho*q1)-mu_eta1/rho+1/2/rho*(sigma_etaa1^2+sigma_etab1^2)...
    -muK/rho+((psiA1*sigmaA)^2+(psiB1*sigmaB)^2)/2/rho...
    -mu_eta1*eta*f(2));



