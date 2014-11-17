function run_TDVP()
%% start ground state calculations
% VMPS_FullSBM(1,0.2,0.1,0)     % VMPS_FullSBM(s,alpha,delta,epsilon)


%% Define GS configuration for TDVP
% Decoupled System-environment:
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141028-1325-SpinBoson-alpha0.05delta0.1epsilon0dk20D5dopt5L49\results.mat');

% iSBM:
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141024-1851-SpinBoson-alpha0.01delta0epsilon0.01dk20D5dopt5L49\results.mat');

% Orth 2010 + Vmat:
% alpha = 0.01
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141025-1342-SpinBoson-alpha0.01delta0.1epsilon0dk20D5dopt5L49\results.mat')
% alpha = 0.05
% load('20141114-2019-SpinBoson-alpha0.05delta0.1epsilon0dk20D5dopt5L49\results.mat')
% alpha = 0.1
% load('20141115-1639-SpinBoson-alpha0.1delta0.1epsilon9dk20D5dopt5L49\results.mat')
% alpha = 0.15
% load('20141115-1639-SpinBoson-alpha0.15delta0.1epsilon9dk20D5dopt5L49\results.mat')
% alpha = 0.2
% load('20141115-1640-SpinBoson-alpha0.2delta0.1epsilon9dk20D5dopt5L49\results.mat')

% Orth 2010 without Vmat:
% alpha = 0.01:
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141025-1750-SpinBoson-alpha0.01delta0.1epsilon0dk20D5dopt5L49\results.mat')
% alpha = 0.05:
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141026-2259-SpinBoson-alpha0.05delta0.1epsilon0dk20D5dopt5L49\results.mat')

% Orthogonal Polynomials, Orth 2010 + Vmat
% alpha = 0.01:
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141114-1625-SpinBoson-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results.mat')
% alpha = 0.05:
% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141114-1617-SpinBoson-alpha0.05delta0.1epsilon0dk20D5dopt5L50\results.mat')
% alpha = 0.2:
load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141116-0229-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L50\results.mat')
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

para.tdvp.filename = sprintf([para.filename(1:end-4),'-Till%dStep%d-noOBBExpand-noBondExpand.mat'],para.tdvp.tmax,para.tdvp.deltaT);

% z-Averaging for log-Discretization
para.tdvp.zAveraging = 0;
if para.tdvp.zAveraging
    para.tdvp.zStep = 0.2;          % Step size for different z values, only for Log-discretization

end

%% Do Time-Evolution with 1-site TDVP
if para.tdvp.zAveraging == 0
    starttime = tic;
    [~, ~, para, results, tmps, tVmat] = tdvp_1site(mps,Vmat,para,results,op);
    results.tdvp.time = toc(starttime);
    save(para.tdvp.filename,'para','Vmat','mps','results','op', 'tmps','tVmat');
    tresults = calTimeObservables(tmps,tVmat,para);
    save(para.tdvp.filename,'para','Vmat','mps','results','op', 'tmps','tVmat','tresults');
else
	basename = para.tdvp.filename;
    for z = 1:-para.tdvp.zStep:1e-5
        %% Set new z-value, prepare chain
        para.z = z;
        [para]=SBM_genpara(para);
        [op,para]=genh1h2term(para);
        [op] = initstorage(mps, Vmat, op,para);
        para.tdvp.filename = sprintf([basename(1:end-4),'-z%.10g.mat'],para.z);

        %% Do time-evolution and save results
        starttime = tic;
        [mps1, Vmat1, para, results, tmps, tVmat] = tdvp_1site(mps,Vmat,para,results,op);
        results.tdvp.time = toc(starttime);
        save(para.tdvp.filename,'para','Vmat','mps','results','op', 'tmps','tVmat');
        tresults = calTimeObservables(tmps,tVmat,para);
        save(para.tdvp.filename,'para','Vmat','mps','results','op', 'tmps','tVmat','tresults');

        %% clear all results from sweep
        clear('tmps','tVmat','mps1','Vmat1');
        results.tdvp = struct();
    end
end
return;

%% calculate observables afterwards:
folders = {'20141115-1639-SpinBoson-alpha0.1delta0.1epsilon9dk20D5dopt5L49',...
           '20141115-1639-SpinBoson-alpha0.15delta0.1epsilon9dk20D5dopt5L49',...
           '20141115-1640-SpinBoson-alpha0.2delta0.1epsilon9dk20D5dopt5L49'};
for k = 1:length(folders)
    files = dir([folders{k},'/','results-Till*']);
    for l = 1:length(files)
        % load file, calculate observables and save it.
        fprintf('Processing %s/%s\n',folders{k},files(l).name);
        load(sprintf('%s/%s',folders{k},files(l).name));
        tresults = calTimeObservables(tmps,tVmat,para);
        save(para.tdvp.filename,'para','Vmat','mps','results','op', 'tmps','tVmat','tresults');
    end
end

end