function filename = VMPS_TreeMPS(varargin)
%% function filename = VMPS_TreeMPS(varargin)
%
%	Initialises the model and starts the TDVP computation with treeMPS.
%	Supports flexible calling syntax via input parser. Input can be Struct with appropriate field names
%
%	Created by FS 26/02/2016
%

%%	Define Input Parser
p = inputParser;					% parses for para
p.CaseSensitive = true;
p.KeepUnmatched = true;
p.FunctionName = 'VMPS_TreeMPS';

pt = inputParser;					% parses for para.tdvp
pt.CaseSensitive = false;
pt.KeepUnmatched = true;
pt.FunctionName = 'VMPS_TreeMPS';

pSBM = inputParser;					% parses for para.SpinBoson
pSBM.CaseSensitive = false;
pSBM.KeepUnmatched = true;
pSBM.FunctionName = 'VMPS_TreeMPS';

%% Define validation functions for string expressions
% para input parser
p.addParameter('model'			,'SpinBoson',@(x) any(validatestring(x,...
								{'SpinBoson','SpinBoson2C', 'DPMES5-7C','testTree'})));		% possible model inputs
p.addParameter('s'				,0	,@isnumeric);
p.addParameter('L'				,10	,@isnumeric);
p.addParameter('alpha'			,0	,@isnumeric);
p.addParameter('delta'			,0	,@isnumeric);
p.addParameter('epsilon'		,0	,@isnumeric);
p.addParameter('dk'				,20	,@isnumeric);
p.addParameter('D'				,5	,@isnumeric);
p.addParameter('d_opt'			,5	,@isnumeric);
p.addParameter('CTshift'		,0	,@isnumeric);
p.addParameter('useDkExpand'	,1	,@isnumeric);
p.addParameter('dkmax'			,1e3,@isnumeric);
p.addParameter('Dmin'			,2	,@isnumeric);
p.addParameter('svmaxtol'	,1e-4	,@isnumeric);
p.addParameter('svmintol'	,10^-4.5,@isnumeric);
p.addParameter('logging'		,0	,@isnumeric);
p.addParameter('rescaling'		,0	,@isnumeric);
p.addParameter('complex'		,0	,@isnumeric);
p.addParameter('tdvp'			,1	,@isnumeric);		% do time-evolution?
p.addParameter('spinposition'	,1	,@isnumeric);
p.addParameter('resume'			,0	,@isnumeric);
p.addParameter('precision'	,5e-15	,@isnumeric);
p.addParameter('parity'			,'n',@isstr);
p.addParameter('eigs_tol'	,1e-8	,@isnumeric);
p.addParameter('foldedChain'	,0	,@isnumeric);
p.addParameter('d_opt_min'		,0	,@isnumeric);
p.addParameter('spinbase'		,'Z',@isstr);
p.addParameter('useshift'		,0	,@isnumeric);
p.addParameter('useTreeMPS'		,1	,@isnumeric);
p.addParameter('useStarMPS'		,0	,@isnumeric);
p.addParameter('useVmat'		,1	,@isnumeric);
p.addParameter('useVtens'		,0	,@isnumeric);

% para.tdvp input parser
% pt.addParameter('model'			,'SpinBoson',@(x) any(validatestring(x,...
% 								{'SpinBoson','SpinBoson2C', 'DPMES5-7C','testTree'})));		% possible model inputs

pt.addParameter('imagT'				,0	,@isnumeric);
pt.addParameter('tmax'				,10	,@isnumeric);
pt.addParameter('deltaT'			,1	,@isnumeric);
pt.addParameter('resume'			,0	,@isnumeric);
pt.addParameter('saveInterval'		,10	,@isnumeric);
pt.addParameter('serialize'			,1	,@isnumeric);
pt.addParameter('logSV'				,0	,@isnumeric);
pt.addParameter('Observables'		,'' ,@isstr);
pt.addParameter('storeMPS'			,0	,@isnumeric);
pt.addParameter('evolveSysTrotter'	,1	,@isnumeric);
pt.addParameter('HEffSplitIsometry'	,1	,@isnumeric);
pt.addParameter('maxExpMDim'		,300,@isnumeric);
pt.addParameter('maxExpVDim'		,700,@isnumeric);
pt.addParameter('expvCustom'		,1	,@isnumeric);
pt.addParameter('useDkExpand'		,1	,@isnumeric);
pt.addParameter('expandOBB'			,1	,@isnumeric);
pt.addParameter('truncateExpandBonds',1	,@isnumeric);
pt.addParameter('maxOBBDim'			,60	,@isnumeric);
pt.addParameter('maxBondDim'		,10	,@isnumeric);
pt.addParameter('zAveraging'		,0	,@isnumeric);
pt.addParameter('logError'			,0	,@isnumeric);
pt.addParameter('logTruncError'		,0	,@isnumeric);
pt.addParameter('rescaling'			,0	,@isnumeric);
% pt.addParameter('storeMPS'			,0	,@isnumeric);

pSBM.addParameter('GroundStateMode'	,'artificial',@(x) any(validatestring(x,...
									{'decoupled','coupled', 'artificial','artTTM'})));			% GroundState preparations
pSBM.addParameter('InitialState'	,'sz',@(x) any(validatestring(x,...
									{'sz', '-sz', 'sx', '-sx', 'sy', '-sy', 'none'})));			% Initial States

%%
p.parse(varargin{:});
para = p.Results;
if para.tdvp == 1
	pt.parse(varargin{:});
	para.tdvp = pt.Results;
	para.complex = 1;
end

%% Defaults to the models

if strfind(para.model,'SpinBoson')
	pSBM.parse(varargin{:});
	para.SpinBoson = pSBM.Results;

%% Setting Hamiltonian for single-spin models
	% H0 = -para.hx./2.*sigmaX -para.hz./2.*sigmaZ
    para.hx = -para.delta;                       % Splitting with sigma_X
    para.hz = -para.epsilon;                     % Splitting with sigma_Z
%% Setting Chain for SBM
	if strcmp(para.model,'SpinBoson')
		para.nChains = 1;
	else
		out = regexp(para.model,'SpinBoson(?<nChains>\d)C','names');
		para.nChains = str2double(out.nChains);
	end

	for mc = 1:para.nChains
		para.chain{mc}.mapping			= 'OrthogonalPolynomials';
		para.chain{mc}.spectralDensity	= 'Leggett_Hard';
		para.chain{mc}.discretization	= 'None';
		% para.chain{1}.discrMethod		= 'Numeric';

		para.chain{mc}.s				= para.s(mc);			% SBM spectral function power law behaviour
		para.chain{mc}.alpha			= para.alpha(mc);		% SBM spectral function magnitude; see Bulla 2003 - 10.1103/PhysRevLett.91.170601
		para.chain{mc}.L				= para.L(mc);
		para.chain{mc}.w_cutoff         = 1;
	end
	
%% Setting TreeMPS for SBM
	para.treeMPS.height    = 1;										% star structure, since only tree node + leaves
	para.treeMPS.maxDegree = para.nChains;
	para.treeMPS.leafIdx   = num2cell((1:para.nChains)'+1);			% maps from chain number to leaf index in para
	para.treeMPS.nodeIdx   = {0+1};									% maps from node number to nodeIdx
%% Set-up parameters for specific ground state preparation!
%     para.SpinBoson.GroundStateMode = 'artificial';				% input arg
        % choose: 'decoupled', 'coupled', 'artificial', 'artTTM'
		% -artificial & artTTM does no optimization! this only sets up an artificial
		%		ground state with <n> = 0 on chain and InitialState 'sz'
		% - 'artificial' for SBM2CT sets up maximally entangled state between odd and even chains in the Bath
		%	odd: thermal bath, even: ancilla
%     para.SpinBoson.InitialState = 'sz';							% input arg
        % choose: 'sz', '-sz', 'sx', '-sx', 'sy', '-sy', 'none'
		% works with all options

	if strcmp(para.SpinBoson.GroundStateMode, 'decoupled')
        para.SpinBoson.t1 = para.t(1);
        para.t(1) = 0;                              % switches off interaction with bath
	end
	if strcmp(para.SpinBoson.InitialState, 'sx') && ~strcmp(para.SpinBoson.GroundStateMode, 'artificial')
        para.SpinBoson.hx = para.hx;
        para.SpinBoson.hz = para.hz;
        para.hz = 0;
        para.hx = 10;                               % some value to do a sx splitting
    elseif strcmp(para.SpinBoson.InitialState, 'sz') && ~strcmp(para.SpinBoson.GroundStateMode, 'artificial')
        para.SpinBoson.hx = para.hx;
        para.SpinBoson.hz = para.hz;
        para.hz = 10;                               % some value to get sz+ eigenstate
        para.hx = 0;
	end
	if isempty(strfind(para.model,'TTM'))
		assert(length(para.spinposition) == 1, 'Only one spin is allowed in SBM');
	end
end



%% Defaults to para
para.loopmax		= 50;						% (minimizeE)
para.increasedk		= 0;						% Tells by how much dk should have been increased to achieve good sv in MPS. start with 0.
para.adjust			= 0;						% (minimizeE) Initialize as 0, no need the edit. To adjust D. Is set = 1 in minimizeE.m
para.dimlock		= 0;						% (minimizeE) set to 0 will change D and d_dop adaptively
para.minDimChange	= 0.01;						% (minimizeE) sets dimlock = 1 if relative dimension change < minDimChange. (larger makes less loops)
para.version		= 'v73';
para.timeslice      = 0;
para.sitej          = 0;

%%%% dk Expansion related %%%%
% only works together with OBB! (para.useVmat = 1)
% para.useDkExpand     = 1;						% input arg. Enable dk expansion, own algorithm. General switch
para.expandprecision = 1e-5;					% unused in dk expand?
para.hasexpanded	 = 0;						% for dk Expand
% Method 1:
para.useDkExpand1	 = 0;						% Expand dk if largest SV of site is smaller than this thershold. EMPIRICAL. Below this value, Vmat seems to need higher dk
para.expandBelowSV   = 0.995;
% Method 2:
para.useDkExpand2	 = 1;						% Expand if wavefunction on site occupies the high energy dimensions
para.dkEx2_tail		 = 0.4;						% tail length [0 1] of occupation in Vmat to analyse
para.dkEx2_maxDev	 = 1.5;						% if std(log10(tail)) < maxDev --> no increase; Measures orders of magnitude in fluctuations of tail.
para.dkEx2_minExp	 = 13;						% if tail below this order than do not expand.

%%%% Boson shift related %%%%
% para.useshift=0;										% input arg
% only choose one of the following Methods
para.useFloShift = 0;                                   % shift all sites if trustsite > 0
para.useFloShift2 = 1;  para.FloShift2minSV = 0.995;    % shift everything if maximum Vmat SV fall below this value
para.useFloShift3 = 0;  para.FloShift3minMaxSV = 1;     % shift every 3rd loop. Shifts only sites where max SV worse than half of the worst SV
para.FloShift3loops = 5;                                % FloShift3 needs logging of results.Vmat_sv!
para.useChengShift=0;                                   % shifts sites < trustsite
para.useEveryShift=0;                                   % shift in every loop

para.shift=zeros(para.nChains,para.L);					% separate shift for each chain!
para.relativeshift=zeros(para.nChains,para.L);
para.relativeshiftprecision=0.01;				% When the relative shift is below this value then stop shifting
% para = maxshift(para);						% TODO!

%% Defaults to para.tdvp
para.tdvp.expvCustomTestAccuracy = 0;		% do expvCustom alongside expv for testing.
para.tdvp.expvCustomTestAccuracyRMS = 0;	% display RMS of expvCustom from expv(); set only if para.tdvp.expvCustomTestAccuracy = 1;
para.tdvp.expvTol = 1e-15;					% error tolerance of expv(); default: 1e-7
para.tdvp.expvM   = 50;						% dim of Krylov subspace in expv(); default: 30
para.tdvp.version = 'v73';

end

function prepareArtState()

end