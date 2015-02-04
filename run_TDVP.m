function run_TDVP(alpha,OBB,Bond)
try
% needed if started from command line:
if isdeployed
	if ischar(alpha), alpha = str2num(alpha); end
	if ischar(OBB), 	OBB 	= str2num(OBB); end
	if ischar(Bond), 	Bond 	= str2num(Bond); end
end

%% start ground state calculations
fileName =  VMPS_FullSBM(1,alpha,0.1,0,200,0);     % VMPS_FullSBM(s,alpha,delta,epsilon,L,rescaling)

%maxNumCompThreads('automatic');			% allows multi-threading in 1pass files
maxNumCompThreads(1);						% safer in terms of results!

load(fileName);

%% Define GS config

% LogZ SBM Orth10 Lambda=2:
% load('20150126-2247-SpinBoson-LogZ-v37TCM66-alpha0.2delta0.1epsilon0dk20D5dopt5L49/results.mat');
% load('20150126-2247-SpinBoson-LogZ-v37TCM66-alpha0.2delta0.1epsilon0dk20D5dopt5L49-artificial/results.mat');

% OrthPol SBM Orth2010 L=50:
% alpha = 0.01
% load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
% load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial/results.mat');

% OrthPol SBM Orth2010 L=200:
% alpha = 0.01
% load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200/results.mat');
% load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial/results.mat');

%% Define GS configuration for TDVP (old)
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

% Orthogonal Polynomials, Orth 2010 + Vmat, L=50 & old
% alpha = 0.01:
% load('20141114-1625-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat')
% load('20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200/results.mat')
% load('20150110-1317-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200/results.mat')
% alpha = 0.05:
% load('20141114-1617-SpinBoson-OrthPol-alpha0.05delta0.1epsilon0dk20D5dopt5L50/results.mat')
% alpha = 0.1:
% load('20141117-0641-SpinBoson-OrthPol-alpha0.1delta0.1epsilon0dk20D5dopt5L50/results.mat')
% alpha = 0.15:
% load('20141117-0642-SpinBoson-OrthPol-alpha0.15delta0.1epsilon0dk20D5dopt5L50/results.mat')
% alpha = 0.2:
% load('20141116-0229-SpinBoson-OrthPol-alpha0.2delta0.1epsilon0dk20D5dopt5L50/results.mat')

% Orthogonal Polynomials, Orth 2010 + OBB, L=200
% if alpha == 0.01
% 	load('20141221-0148-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200/results.mat')
% elseif alpha == 0.05
% 	load('20141221-0148-SpinBoson-OrthPol-alpha0.05delta0.1epsilon0dk20D5dopt5L200/results.mat')
% elseif alpha == 0.1
% 	load('20141221-0151-SpinBoson-OrthPol-alpha0.1delta0.1epsilon0dk20D5dopt5L200/results.mat')
% elseif alpha == 0.15
% 	load('20141221-0151-SpinBoson-OrthPol-alpha0.15delta0.1epsilon0dk20D5dopt5L200/results.mat')
% elseif alpha == 0.2
% 	load('20141221-0153-SpinBoson-OrthPol-alpha0.2delta0.1epsilon0dk20D5dopt5L200/results.mat')
% end

%% Define TDVP parameters
para.tdvp.tmax = 1e5;
para.tdvp.tmin = 0.1;
    % For PPC:
    %   H defined in eV, h\bar left out
    %   -> real tmax = T * 6.58211928(15)×10^-16
para.tdvp.deltaT = 0.1;                 % size of timeslice in units:
para.tdvp.timescale = 'Exp10';			% Exponential time steps
% para.tdvp.t = 0:para.tdvp.deltaT:para.tdvp.tmax;
para.tdvp.t = [0,10.^(log10(para.tdvp.tmin):para.tdvp.deltaT:log10(para.tdvp.tmax))];		% for use with log10
para.tdvp.maxExpMDim = 260;			% For Lappy: 100, OE-PC: 80, pc52: 260; System dependent, use benchmark!
para.tdvp.maxExpVDim = 800;				% higher dim -> use expvCustom() if expvCustom == 1. Number from benchmarking. Lappy: 600, Haswell: 800; maxExpMDim < maxExpVDim
para.tdvp.expvCustom = 1;				% 1 for Custom programmed, 0 for standard expv()
para.tdvp.expvCustomTestAccuracy = 0;	% do expvCustom alongside expv for testing.
para.tdvp.expvCustomTestAccuracyRMS = 0;	% display RMS of expvCustom from expv(); set only if para.tdvp.expvCustomTestAccuracy = 1;
para.tdvp.expvTol = 1e-15;              % error tolerance of expv(); default: 1e-7
para.tdvp.expvM   = 50;                 % dim of Krylov subspace in expv(); default: 30
    % Sets threshold size for matrix exponential:
    %   if dim(A) < : use built-in expm(At)*v
    %   else        : use Expokit expv(t,A,v, expvTol, expvM)
    %   set maxExpMDim = 0 to only use expv()
% OBB settings
para.tdvp.expandOBB = OBB;
% Bond-Dim settings
para.tdvp.truncateExpandBonds = Bond;
% Calculate max Bond Dim: 1GB for array (l,r,n,l,r,n) with n around 20,
% 1 complex double needs 16byte. -> 20^6 * 16byte < 1GB
para.tdvp.maxBondDim = 15;
% z-Averaging for log-Discretization
para.tdvp.zAveraging = 0;
if para.tdvp.zAveraging
    para.tdvp.zStep = 0.2;          % Step size for different z values, only for Log-discretization
end

%% Format Filename
para.tdvp.version = 'v39';
if isfield(para.tdvp,'filename')
	%% Continued TDVP
	para.tdvp.fromFilename = para.tdvp.filename;		% save the reference to continued file
	% input the TDVP matrices!
	mps = tmps;
	Vmat = tVmat;
else
	para.tdvp.fromFilename = para.filename;				% or reference to VMPS Ground State File!
end

if strcmp(para.tdvp.timescale,'Exp10')
	para.tdvp.filename = sprintf([para.filename(1:end-4),'-Till%dExpStep%.2g%s.mat'],para.tdvp.tmax,para.tdvp.deltaT,para.tdvp.version);
else
	para.tdvp.filename = sprintf([para.filename(1:end-4),'-Till%dStep%.2g%s.mat'],para.tdvp.tmax,para.tdvp.deltaT,para.tdvp.version);
end
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
para.tdvp.filename = sprintf([para.tdvp.filename(1:end-4),'-%dcore.mat'],maxNumCompThreads);

%% Check for valid simulation arguments and declare unsettables
if para.tdvp.expvCustomTestAccuracyRMS
	assert(para.tdvp.expvCustomTestAccuracy == 1, 'Only set para.tdvp.expvCustomTestAccuracyRMS == 1 if para.tdvp.expvCustomTestAccuracy == 1');
end
if para.tdvp.zAveraging
	assert(strcmp(para.chainMapping,'LogDiscrZitko'),'Use z-Averaging only with Logarithmic Discretization!');
end
if para.tdvp.expvCustom
	assert(para.tdvp.maxExpMDim <= para.tdvp.maxExpVDim,'maxExpMDim <= maxExpVDim ! Everything else has no sense.');
end

para.tdvp.expvCustomNow = 0;			% only used inside the program
para.tdvp.rescaling = 0;                % turn on/off rescaling in TDVP; needs to be off since H in exponent!
para.rescaling = para.tdvp.rescaling;
para.complex = 1;						% necessary since time-evolution is complex

%% Copy to scratch for computation
if ~strcmp(computer,'PCWIN64')
	[~, name] = system('hostname');
	para.tdvp.hostname = strtrim(name);
	save(sprintf([para.tdvp.filename(1:end-4),'-incomplete-%s.mat'],para.tdvp.hostname),'para','results');
	tempDir = '/scratch/fayns2/TDVPtemp/';tempFold = fileparts(para.filename);
	currentDir = pwd;
	addpath(currentDir);
	mkdir([tempDir,tempFold]);
	copyfile(para.filename,[tempDir,tempFold]);
	cd(tempDir);
end

%% Do Time-Evolution with 1-site TDVP
if para.tdvp.zAveraging == 0
    para.tdvp.starttime = tic;
    tdvp_1site(mps,Vmat,para,results,op);
    load(para.tdvp.filename,'para','tmps','tVmat','results');
	results.tdvp.time = toc(para.tdvp.starttime);
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
        para.tdvp.starttime = tic;
        tdvp_1site(mps,Vmat,para,results,op);
        load(para.tdvp.filename,'para','tmps','tVmat','results');
		results.tdvp.time = toc(para.tdvp.starttime);
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
	%delete([currentDir,'/',para.tdvp.filename(1:end-4),'-incomplete.mat']);
	sendmailCAM('fayns2@cam.ac.uk',...
         'TDVP job completed',sprintf('The job \n %s\nHas successfully completed.',para.tdvp.filename));
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
catch err
	fprintf([getReport(err),'\n']);
	sendmailCAM('fayns2@cam.ac.uk',...
         'TDVP job ERROR',sprintf('The job \n %s\nHas encountered an error:\n%s',para.tdvp.filename,getReport(err)));
end
end