function run_TDVP(s,alpha,OBB,Bond)
% try
% needed if started from command line:
if isdeployed
	if ischar(s), s = str2num(s); end
	if ischar(alpha), alpha = str2num(alpha); end
	if ischar(OBB), 	OBB 	= str2num(OBB); end
	if ischar(Bond), 	Bond 	= str2num(Bond); end
end

%% start ground state calculations
fileName =  VMPS_FullSBM(s,alpha,0.1,0,100,0);     % VMPS_FullSBM(s,alpha,delta,epsilon,L,rescaling)

%maxNumCompThreads('automatic');			% allows multi-threading in 1pass files
maxNumCompThreads(1);						% safer in terms of results!

%% Define GS config

load(fileName);

% load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150327-1434-SpinBoson-OrthPol-v42TCMde10-s0.5-alpha0.01delta0.1epsilon0dk20D5dopt5L100\results.mat')
% load(sprintf('20150512-1604-SpinBoson-OrthPol-v43TCMde10-alpha%gdelta0.1epsilon0dk30D5dopt5L300-artificial/results-Till500Step0.2v43-OBBExpand-noBondExpand-expvCustom800-1core-small.mat',alpha));
% load('20150723-1320-SpinBoson2C-OrthPol-v49-alpha0.2delta0.1epsilon0dk6D5dopt5L10/results.mat');
% load('20150723-1501-SpinBoson2folded-OrthPol-v49-alpha0.2delta0.1epsilon0dk36D5dopt5L10/results.mat');

% Kast 2013 Fig 4, s=0.75 OrthPol
% load(sprintf('20150307-0341-SpinBoson-OrthPol-v41TCMde9-s0.75-alpha%gdelta0.1epsilon0dk20D5dopt5L50-artificial/results.mat',alpha));

%% Define TDVP parameters
para.tdvp.tmax = 240;
para.tdvp.deltaT = 0.2;                % size of timeslice in units:
    % For PPC:
    %   H defined in eV, h\bar left out
    %   -> real tmax = T * 6.58211928(15)×10^-16
para.tdvp.t = 0:para.tdvp.deltaT:para.tdvp.tmax;
para.tdvp.resume = 0;					% additionally control if want to resume!
para.tdvp.saveInterval = 5;				% save '-small.mat' every n-th step
para.tdvp.serialize = 0;				% much faster I/O saving
para.tdvp.logSV = 0;					% if 1 only log SV, if 0 only log vNE (saves mem) if -1 log none!
% para.tdvp.extractStarInterval = 1;		% in [t]; for calculating star occupation! Comment if not needed!
para.tdvp.extractObsInterval  = 0.2;		% in [t]; mod(extractStarInterval, extractObsInterval) = 0 !! extractObsInterval = n*deltaT
para.tdvp.Observables = '.n.j.s.sn.sx.';% n: occupation, j: current, s: spin, sn: star n, sx: star polaron
para.tdvp.storeMPS = 0;					% save tmps or not!
para.tdvp.maxExpMDim = 100;				% For Lappy: 100, OE-PC: 80, pc52: 260; System dependent, use benchmark!
para.tdvp.maxExpVDim = 600;				% higher dim -> use expvCustom() if expvCustom == 1. Number from benchmarking. Lappy: 600, Haswell: 800; E5: 960 maxExpMDim < maxExpVDim
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
para.tdvp.maxBondDim = 20;
para.tdvp.maxOBBDim  = 20;
para.svmaxtol = 10^-4;
para.svmintol = 10^-4.5;
% z-Averaging for log-Discretization
para.tdvp.zAveraging = 0;
if para.tdvp.zAveraging
    para.tdvp.zStep = 0.2;          % Step size for different z values, only for Log-discretization
	if ~exist('stepFrom','var')
		stepFrom = 1;					% comment if need override!
	end
end

tresults = [];						% empty variable initialization
%% Format Filename
para.tdvp.version = 'v50';
if isfield(para.tdvp,'filename')
	%% Continued TDVP remember filename to load after directory change!
	% from File can be -small.mat!
	para.tdvp.fromFilename = para.tdvp.filename;		% save the reference to continued file
	% input the TDVP matrices!
% 	mps = tmps;
% 	Vmat = tVmat;
else
	para.tdvp.fromFilename = para.filename;				% or reference to VMPS Ground State File!
end

para.tdvp.filename = sprintf([para.filename(1:end-4),'-Till%dStep%.2g%s.mat'],para.tdvp.tmax,para.tdvp.deltaT,para.tdvp.version);

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
para.tdvp.filenameSmall = [para.tdvp.filename(1:end-4),'-small.mat'];		% only for para, tresults
% Set MPS filename if needed
if para.tdvp.storeMPS == 1
	para.tdvp.filenameMPS = [para.tdvp.filename(1:end-4),'-MPS.mat'];		% only for tmps, tVmat
end

%% Check for valid simulation arguments and declare unsettables
if para.tdvp.expvCustomTestAccuracyRMS
	assert(para.tdvp.expvCustomTestAccuracy == 1, 'Only set para.tdvp.expvCustomTestAccuracyRMS == 1 if para.tdvp.expvCustomTestAccuracy == 1');
end
if para.tdvp.zAveraging
	assert(strcmp(para.chain.discretization,'LogZ'),'Use z-Averaging only with Logarithmic Discretization!');
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
	para.tdvp.scratchDir = '/scratch/fayns2/TDVPtemp/'; tempFold = fileparts(para.filename);
	currentDir = pwd;
	addpath(currentDir);
	if ~exist(para.tdvp.scratchDir,'dir')
		mkdir(para.tdvp.scratchDir);
	end
	mkdir([para.tdvp.scratchDir,tempFold]);
	copyfile(para.filename,[para.tdvp.scratchDir,tempFold]);
	save(sprintf([para.tdvp.filename(1:end-4),'-incomplete-%s.mat'],para.tdvp.hostname),'para','results');
	cd(para.tdvp.scratchDir);
end

if exist(para.tdvp.filename,'file') && para.tdvp.resume
	% either load and resume or ignore?
	try
		load(para.tdvp.filename);	% all the rest
	catch
		disp('MAT-file corrupt, load BAK');
		load([para.tdvp.filename(1:end-4),'.bak'],'-mat');		% use backup if file corrupted
	end
	if isnumeric(para)
		% everything was serialized!
		para     = hlp_deserialize(para);
		op       = hlp_deserialize(op);
		results  = hlp_deserialize(results);
		tresults = hlp_deserialize(tresults);
		mps      = hlp_deserialize(mps);
		Vmat     = hlp_eserialize(Vmat);
	end
	para.tdvp.resume = 1;		% need override to restore value!
elseif ~exist(para.tdvp.filename,'file') && isempty(strfind(para.tdvp.fromFilename,'results.mat'))
	% Only load fromFilename if it was a TDVP file!
	% then also para and tresults are already loaded!
	% important: load mps and Vmat!
	load(para.tdvp.fromFilename, 'results','op','mps','Vmat');
	if isnumeric(results)
		% everything was serialized!
		op       = hlp_deserialize(op);
		results  = hlp_deserialize(results);
		mps      = hlp_deserialize(mps);
		Vmat     = hlp_eserialize(Vmat);
	end
end

%% Do Time-Evolution with 1-site TDVP
if para.tdvp.zAveraging == 0
    para.tdvp.starttime = tic;
    tdvp_1site(mps,Vmat,para,results,op,tresults);
else
	basename = para.tdvp.filename;
    for z = stepFrom:-para.tdvp.zStep:1e-5
        %% Set new z-value, prepare chain
        para.z = z;
        [para]=SBM_genpara(para);
        [op,para]=genh1h2term(para);
        [op] = initstorage(mps, Vmat, op,para);
        para.tdvp.filename = sprintf([basename(1:end-4),'-z%.10g.mat'],para.z);

        %% Do time-evolution and save results
        para.tdvp.starttime = tic;
        tdvp_1site(mps,Vmat,para,results,op);

		if ~strcmp(computer,'PCWIN64')
			copyfile(para.tdvp.filenameSmall,[currentDir,'/',para.tdvp.filenameSmall]);
		end
        %% clear all results from sweep
        clear('tmps','tVmat','mps1','Vmat1');
        results.tdvp = struct();
    end
end
if ~strcmp(computer,'PCWIN64')
	copyfile(para.tdvp.filenameSmall,[currentDir,'/',para.tdvp.filenameSmall]);
	%delete([currentDir,'/',para.tdvp.filename(1:end-4),'-incomplete.mat']);
	sendmailCAM('fayns2@cam.ac.uk',...
         'TDVP job completed',sprintf('The job \n %s\nHas successfully completed.',para.tdvp.filename));
	exit;
end

% catch err
% 	fprintf([getReport(err),'\n']);
% 	sendmailCAM('fayns2@cam.ac.uk',...
%          'TDVP job ERROR',sprintf('The job \n %s\nHas encountered an error:\n%s',para.tdvp.filename,getReport(err)));
% end
end