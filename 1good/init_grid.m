function grid = init_grid(s)
% Initializes grid to hold data about the equilibrium, state space, etc.
%
% Written by William Chen, Mar. 2019

if not(isfield(s, 'start'))
    s.start = 1e-3; % perturb from zero b/c singularity
end

% make eta state space
if isfield(s, 'etagrid') == 1
    grid.eta = s.etagrid; % pre-specified eta vector
elseif isfield(s, 'chebnodes') == 1
    if s.chebnodes == 1
        grid.eta = makegrid(s.start, s.end - s.start, s.N);
    end
else % construct the state space automatically
    if isfield(s, 'N')
        if not(s.N / 2 == ceil(s.N/2))
            error('Size of state space is not even.')
        end
        grid.eta = linspace(s.start, .5, s.N/2)';
    else
        grid.eta = linspace(s.start, .5, 125)';
    end
    % apply twisting to get uneven state space
    if not(isfield(s, 'twist'))
        s.twist = 2;
    end
    grid.eta(1:s.N/2) = grid.eta(1:s.N/2).^s.twist / .5.^(s.twist - 1);
    grid.eta(s.N/2)   = (grid.eta(s.N/2) - grid.eta(s.N/2-1)) * 2/3 + grid.eta(s.N/2-1);
    grid.eta          = [grid.eta; flipud(1 - grid.eta(1:s.N/2))]; % symmetric grid
end
grid.dim1 = length(grid.eta);

% Make vectors for holding data
grid.mu_eta = zeros(grid.dim1, 1);
grid.mu_vAerr = zeros(grid.dim1, 1);
grid.mu_vBerr = zeros(grid.dim1, 1);
grid.sig_etaA = zeros(grid.dim1, 1);
grid.sig_etaB = zeros(grid.dim1, 1);
grid.sig_qA = zeros(grid.dim1, 1);
grid.sig_qB = zeros(grid.dim1, 1);
grid.sig_xiA = zeros(grid.dim1, 1);
grid.sig_xiB = zeros(grid.dim1, 1);
grid.sig_zetaA = zeros(grid.dim1, 1);
grid.sig_zetaB = zeros(grid.dim1, 1);
grid.sig_vAA = zeros(grid.dim1, 1);
grid.sig_vAB = zeros(grid.dim1, 1);
grid.sig_vBA = zeros(grid.dim1, 1);
grid.sig_vBB = zeros(grid.dim1, 1);
grid.psiAa = zeros(grid.dim1, 1);
grid.psiAb = zeros(grid.dim1, 1);
grid.psiBa = zeros(grid.dim1, 1);
grid.psiBb = zeros(grid.dim1, 1);
grid.psiA = zeros(grid.dim1, 1);
grid.psiB = zeros(grid.dim1, 1);
grid.q = zeros(grid.dim1, 1);
grid.qp = zeros(grid.dim1, 1);
grid.qpp = zeros(grid.dim1, 1);
grid.mu_q = zeros(grid.dim1, 1);
grid.muK = zeros(grid.dim1, 1);
grid.Pa = zeros(grid.dim1, 1);
grid.Pb = zeros(grid.dim1, 1);
end