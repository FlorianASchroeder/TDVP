function run_TDVP(alpha)
%% start ground state calculations
fileName =  VMPS_FullSBM(1,alpha,0.1,0)     % VMPS_FullSBM(s,alpha,delta,epsilon)

% load(fileName);

%% Define GS configuration for TDVP
% Decoupled System-environment:
% load('20141028-1325-SpinBoson-alpha0.05delta0.1epsilon0dk20D5dopt5L49/results.mat');

% iSBM:
% load('20141024-1851-SpinBoson-alpha0.01delta0epsilon0.01dk20D5dopt5L49/results.mat');

% Orth 2010 + Vmat:
% alpha = 0.01
% load('20141025-1342-SpinBoson-alpha0.01delta0.1epsilon0dk20D5dopt5L49/results.mat')
% alpha = 0.05
% load('20141114-2019-SpinBoson-alpha0.05delta0.1epsilon0dk20D5dopt5L49\results.mat')
% alpha = 0.1
% load('20141117-0405-SpinBoson-alpha0.1delta0.1epsilon0dk20D5dopt5L49/results.mat')
% alpha = 0.15
% load('20141117-0406-SpinBoson-alpha0.15delta0.1epsilon0dk20D5dopt5L49/results.mat')
% alpha = 0.2
% load('20141117-0406-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L49/results.mat')
% load('20141117-0531-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L84/results.mat')

% Orth 2010 without Vmat:
% alpha = 0.01:
% load('20141025-1750-SpinBoson-alpha0.01delta0.1epsilon0dk20D5dopt5L49/results.mat')
% alpha = 0.05:
% load('20141026-2259-SpinBoson-alpha0.05delta0.1epsilon0dk20D5dopt5L49/results.mat')

% Orthogonal Polynomials, Orth 2010 + Vmat
% alpha = 0.01:
load('20141114-1625-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat')
% load('20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200/results.mat')

% alpha = 0.05:
% load('20141114-1617-SpinBoson-OrthPol-alpha0.05delta0.1epsilon0dk20D5dopt5L50/results.mat')
% alpha = 0.1:
% load('20141117-0641-SpinBoson-OrthPol-alpha0.1delta0.1epsilon0dk20D5dopt5L50/results.mat')
% alpha = 0.15:
% load('20141117-0642-SpinBoson-OrthPol-alpha0.15delta0.1epsilon0dk20D5dopt5L50/results.mat')
% alpha = 0.2:
% load('20141116-0229-SpinBoson-OrthPol-alpha0.2delta0.1epsilon0dk20D5dopt5L50/results.mat')

%% Define TDVP parameters
para.tdvp.tmax = 20;
    % For PPC:
    %   H defined in eV, h\bar left out
    %   -> real tmax = T * 6.58211928(15)�10^-16
para.tdvp.deltaT = 4;                 % size of timeslice in units:
para.tdvp.t = 0:para.tdvp.deltaT:para.tdvp.tmax;
para.tdvp.maxExpMDim = 10^2;			% For Lappy: 100, OE-PC: 80, pc52: 260; System dependent, use benchmark!
para.tdvp.maxExpVDim = 800;				% higher dim -> use expvCustom() if expvCustom == 1. Number from benchmarking. Lappy: 600, Haswell: 800
para.tdvp.expvCustom = 1;				% 1 for Custom programmed, 0 for standard expv()
para.tdvp.expvCustomNow = 0;			% only set inside the program
para.tdvp.expvCustomTestAccuracy = 0;	% do expvCustom alongside expv for testing.
para.complex = 1;
para.tdvp.expvTol = 1e-15;              % error tolerance of expv(); default: 1e-7
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

%% Format Filename
if isfield(para.tdvp,'filename')
	%% Continued TDVP
	para.tdvp.fromFilename = para.tdvp.filename;		% save the reference to continued file
	% input the TDVP matrices!
	mps = tmps;
	Vmat = tVmat;
end

para.tdvp.filename = sprintf([para.filename(1:end-4),'-Till%dStep%d.mat'],para.tdvp.tmax,para.tdvp.deltaT);
if para.tdvp.expandOBB
	para.tdvp.filename = sprintf([para.tdvp.filename(1:end-4),'-OBBExpand.mat']);
else
	para.tdvp.filename = sprintf([para.tdvp.filename(1:end-4),'-noOBBExpand.mat']);
end
if para.tdvp.truncateExpandBonds
	para.tdvp.filename = sprintf([para.tdvp.filename(1:end-4),'-BondExpand%d.mat'],para.tdvp.maxBondDim);
else
	para.tdvp.filename = sprintf([para.tdvp.filename(1:end-4),'-noBondExpand.mat']);
end
if para.tdvp.expvCustom
	para.tdvp.filename = sprintf([para.tdvp.filename(1:end-4),'-expvCustom%d.mat'],para.tdvp.maxExpVDim);
end

% z-Averaging for log-Discretization
para.tdvp.zAveraging = 0;
if para.tdvp.zAveraging
    para.tdvp.zStep = 0.2;          % Step size for different z values, only for Log-discretization

end

%% Copy to scratch for computation
if ~strcmp(computer,'PCWIN64')
	save([para.tdvp.filname(1:end-4),'-incomplete.mat'],'para','results');
	tempDir = '/scratch/fayns2/TDVPtemp/';tempFold = fileparts(para.filename);
	currentDir = pwd;
	addpath(currentDir);
	mkdir([tempDir,tempFold]);
	copyfile(para.filename,[tempDir,tempFold]);
	cd(tempDir);
end

%% Do Time-Evolution with 1-site TDVP
if para.tdvp.zAveraging == 0
    starttime = tic;
    tdvp_1site(mps,Vmat,para,results,op);
    load(para.tdvp.filename,'para','tmps','tVmat','results');
	results.tdvp.time = toc(starttime);
    tresults = calTimeObservables(tmps,tVmat,para);
    save(para.tdvp.filename,'tresults','-append');
	save([para.tdvp.filename(1:end-4),'-small.mat'],'para','results','tresults');
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
        tdvp_1site(mps,Vmat,para,results,op);
        load(para.tdvp.filename,'para','tmps','tVmat','results');
		results.tdvp.time = toc(starttime);
        tresults = calTimeObservables(tmps,tVmat,para);
        save(para.tdvp.filename,'tresults','-append');
		save([para.tdvp.filename(1:end-4),'-small.mat'],'para','results','tresults');

		if ~strcmp(computer,'PCWIN64')
			copyfile([para.tdvp.filename(1:end-4),'-small.mat'],[currentDir,'/',para.tdvp.filename(1:end-4),'-small.mat']);
		end
        %% clear all results from sweep
        clear('tmps','tVmat','mps1','Vmat1');
        results.tdvp = struct();
    end
end
if ~strcmp(computer,'PCWIN64')
	copyfile([para.tdvp.filename(1:end-4),'-small.mat'],[currentDir,'/',para.tdvp.filename(1:end-4),'-small.mat']);
	exit;
end
return;

%% calculate observables afterwards in z:
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