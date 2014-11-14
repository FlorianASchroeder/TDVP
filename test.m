%% start ground state calculations
% VMPS_FullSBM(1,0.01,0.1,0)     % VMPS_FullSBM(s,alpha,delta,epsilon)


%% Define GS configuration for TDVP
% Decoupled System-environment:
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141028-1325-SpinBoson-alpha0.05delta0.1epsilon0dk20D5dopt5L49\results.mat');

% iSBM:
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141024-1851-SpinBoson-alpha0.01delta0epsilon0.01dk20D5dopt5L49\results.mat');

% Orth 2010 + Vmat:
% alpha = 0.01
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141025-1342-SpinBoson-alpha0.01delta0.1epsilon0dk20D5dopt5L49\results.mat')

% Orth 2010 without Vmat:
% alpha = 0.01:
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141025-1750-SpinBoson-alpha0.01delta0.1epsilon0dk20D5dopt5L49\results.mat')
% alpha = 0.05:
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141026-2259-SpinBoson-alpha0.05delta0.1epsilon0dk20D5dopt5L49\results.mat')

% Orthogonal Polynomials, Orth 2010 + Vmat
% alpha = 0.01:
load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141114-1625-SpinBoson-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results.mat')
% alpha = 0.05:
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141114-1617-SpinBoson-alpha0.05delta0.1epsilon0dk20D5dopt5L50\results.mat')
%% Define TDVP parameters
para.tdvp.tmax = 325;
    % For PPC:
    %   H defined in eV, h\bar left out
    %   -> real tmax = T * 6.58211928(15)×10^-16
para.tdvp.deltaT = 4;                 % size of timeslice in units:
para.tdvp.t = 0:para.tdvp.deltaT:para.tdvp.tmax;
para.tdvp.maxExpMDim = 10^2;            % For Lappy: 100, OE-PC: 80, pc52: 260; System dependent, use benchmark!
para.tdvp.expvTol = 1e-15;               % error tolerance of expv(); default: 1e-7
para.tdvp.expvM   = 50;                 % dim of Krylov subspace in expv(); default: 30
    % Sets threshold size for matrix exponential:
    %   if dim(A) < : use built-in expm(At)*v
    %   else        : use Expokit expv(t,A,v, expvTol, expvM)
    %   set maxExpMDim = 0 to only use expv()
para.tdvp.rescaling = 0;                % turn on/off rescaling in TDVP
para.rescaling = para.tdvp.rescaling;
% OBB settings
para.tdvp.expandOBB = 0;
% Bond-Dim settings
para.tdvp.truncateExpandBonds = 0;
% Calculate max Bond Dim: 1GB for array (l,r,n,l,r,n) with n around 20,
% 1 complex double needs 16byte. -> 20^6 * 16byte < 1GB
para.tdvp.maxBondDim = 20;

para.tdvp.filename = sprintf([para.filename(1:end-4),'-Till%dStep%d-noOBBorBondExpand.mat'],para.tdvp.tmax,para.tdvp.deltaT);

% z-Averaging for log-Discretization
para.tdvp.zAveraging = 0;
if para.tdvp.zAveraging
    para.tdvp.zStep = 0.2;          % Step size for different z values, only for Log-discretization

end

%% Do Time-Evolution with 1-site TDVP
if para.tdvp.zAveraging == 0
    [mps1, Vmat1, para, results, tmps, tVmat] = tdvp_1site(mps,Vmat,para,results,op);

    save(para.tdvp.filename,'para','Vmat','mps','results','op', 'tmps','tVmat');
else
    for z = 0.9:-para.tdvp.zStep:1e-5
        %% Set new z-value, prepare chain
        para.z = z;
        [para]=SBM_genpara(para);
        [op,para]=genh1h2term(para);
        [op] = initstorage(mps, Vmat, op,para);
        para.tdvp.filename = sprintf([para.tdvp.filename(1:end-4),'z%.10g-OBBExpand-noBondExpand.mat'],para.z);

        %% Do time-evolution and save results
        [mps1, Vmat1, para, results, tmps, tVmat] = tdvp_1site(mps,Vmat,para,results,op);
        save(para.tdvp.filename,'para','Vmat','mps','results','op', 'tmps','tVmat');

        %% clear all results from sweep
        clear('tmps','tVmat','mps1','Vmat1');
        results.tdvp = struct();
        para.tdvp.filename = sprintf([para.filename(1:end-4),'-Till%dStep%d.mat'],para.tdvp.tmax,para.tdvp.deltaT);

    end
end

%% calculate observables:
fprintf('Calculate Observables:\n');

tresults.nx = zeros(size(tmps,1),para.L);

for i = 1:size(tmps,1)
    fprintf('%g-',i);

    % 1. Chain Occupation
    tresults.nx(i,:) = getObservable({'occupation'},tmps(i,:),tVmat(i,:),para);

    if strcmp(para.model, 'SpinBoson')
        % 1. Spin Observables
        if ~isfield(tresults,'spin')
            tresults.spin.sx = zeros(size(tmps,1),1);
            tresults.spin.sy = zeros(size(tmps,1),1);
            tresults.spin.sz = zeros(size(tmps,1),1);
            tresults.spin.visibility = zeros(size(tmps,1),1);
        end
        temp = getObservable({'spin'},tmps(i,:),tVmat(i,:),para);
        tresults.spin.sx(i) = temp.sx;
        tresults.spin.sy(i) = temp.sy;
        tresults.spin.sz(i) = temp.sz;
        tresults.spin.visibility(i) = sqrt(temp.sx^2+temp.sy^2);
    end


    if strcmp(para.model, 'MLSBM')
        % 2. PPC Wavefunction
        if ~isfield(tresults,'PPCWavefunction')
            tresults.PPCWavefunction = zeros(size(tmps,1),16);
        end
        tresults.PPCWavefunction(i,:) = diag(getObservable({'rdm',1},tmps(i,:),tVmat(i,:),para));

        % 3. Participation on ring
        if ~isfield(tresults,'participation')
            tresults.participation = zeros(size(tmps,1),1);
        end
        tresults.participation(i) = getObservable({'participation'},tmps(i,:),tVmat(i,:),para);
    end
end
save(para.tdvp.filename,'para','Vmat','mps','results','op', 'tmps','tVmat','tresults');
fprintf('\n');
