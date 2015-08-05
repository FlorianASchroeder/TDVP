function fileName = VMPS_FullSBM(s,alpha,delta,epsilon,L,rescaling)
%Variational matrix product method to study spin-boson model. Rok's
%logrithimic discretization and optimal boson basis are implemented in
%this code.
% adapted to use the SBM as benchmark system

% modified to support a multi-level spin site with choosable coupling to the bath.


%Cheng Guo @Munich
% 3 Sep 2010


% mod by Florian Schroeder @Cambridge
% 23 October 2014 -
% copied from MLSBM 23/10/14


maxNumCompThreads(1);
format short e

starttime = tic;
if isdeployed           % take care of command line arguments
% 	if ischar(hx), hx = str2num(hx); end
% 	if ischar(hz), hz = str2num(hz); end
	if ischar(s), s = str2num(s); end
	if ischar(alpha), alpha = str2num(alpha); end
	if ischar(delta), delta = str2num(delta); end
	if ischar(epsilon), epsilon = str2num(epsilon); end
	if ischar(L), L = str2num(L); end
	if ischar(rescaling), rescaling = str2num(rescaling); end
% 	if ischar(Lambda), Lambda = str2num(Lambda); end
% 	if ischar(parity), parity = str2num(parity); end
end

%% Choose model and chain mapping
para.model='SpinBoson';
    % choose: 'SpinBoson', 'SpinBoson2folded', 'MLSpinBoson','ImpurityQTN'
	%         '2SpinPhononModel', 'SpinBoson2C'
% para.chainMapping = 'OrthogonalPolynomials';
para.nEnvironments   = 1;
	% number of different spectral functions
	% supported 1 to any
para.nChains		 = 1;
	% number of chains
	% 1 for folded, can have nEnvironments = 2;
	% = nEnvironments for multi-chain models;
%% System Definitions:
if ~strcmp(para.model,'MLSpinBoson') && ~strcmp(para.model,'2SpinPhononModel')
	% setting para for single-spin models
    para.hx = -delta;                       % Splitting with sigma_X
    para.hz = -epsilon;                     % Splitting with sigma_Z
end

%% Chain Definitions:
% para.chain{i}.mapping
	% choose: 'OrthogonalPolynomials','LanzcosTriDiag', Stieltjes'
	%			- 'LanzcosTriDiag' Lanczos tridiagonalization
	%			- 'Stieltjes': needs discretization, based on Orthogonal
	%						   Polynomials, less stable but faster than Lanczos
	%			- 'OrthogonalPolynomials'
% para.chain{i}.spectralDensity
	% choose: 'Leggett_Hard', 'Leggett_Soft', 'Renger'
	% Leggett for SBM;  see DOI: 10.1103/PhysRevLett.91.170601
	% Renger for MLSBM; see DOI: 10.1016/j.bbabio.2012.02.016
% para.chain{i}.discretization
	% choose: 'None', 'Linear', 'LogZ'
	%	'None': exact mapping
	%	other methods not yet supported
para.L = L;
%% chain 1:
para.chain{1}.mapping			= 'OrthogonalPolynomials';
para.chain{1}.spectralDensity	= 'Leggett_Hard';
para.chain{1}.discretization	= 'None';
	% para.chain{1}.discrMethod		= 'Numeric';

para.chain{1}.s					= s;			% SBM spectral function power law behaviour
para.chain{1}.alpha				= alpha;		% SBM spectral function magnitude; see Bulla 2003 - 10.1103/PhysRevLett.91.170601
para.chain{1}.L					= para.L;
if alpha == 0 && para.chain{1}.L == 0
	para.chain{1}.L = 10;						% otherwise encounter error
end

%% chain 2:
% para.chain{2}					= para.chain{1};		% simple copy
% para.chain{2}.mapping			= 'OrthogonalPolynomials';
% para.chain{2}.spectralDensity	= 'Leggett_Hard';

% para.chain{3} = para.chain{2};
assert(para.nEnvironments == length(para.chain),'number of environments is wrong');		% redundant, sanity check!
%% Parameters
for k = 1:para.nEnvironments
	if strcmp(para.chain{k}.mapping,'OrthogonalPolynomials')
		%% now only for SBM, to be extended for any J(w)
		% no need to define:
		% z, Lambda
		% since site energies converge to w_c/2 = 0.5, Optimum chain length can
		% not easily be determined -> give para.L
		% based on analytic results only! For numeric approach use Stieltjes

		if para.chain{k}.L == 0					% chain length override
			para.chain{k}.L = 200;				% default chain length if input L=0
		end
		para.rescaling = 0;						% only for LogZ discretization applicable
	% 	para.Lambda = 2;						% Bath log Discretization parameters in case rescaling = 1
	%   para.z	    = 1;
	% 	if para.Lambda == 1
	% 		assert(para.rescaling == 0, 'Please switch off rescaling when using Lambda = 1');
	% 		assert(~strcmp(para.chain.discretization,'LogZ'), 'Lambda = 1 not possible with LogZ discretization!');
	% 	end

	elseif strcmp(para.chain{k}.mapping,'LanzcosTriDiag') || strcmp(para.chain{k}.mapping, 'Stieltjes')
		% 	if Lambda > 1 -> LogZ
		% 	if Lambda = 1 -> Linear discretization, rescaling = 0
		para.chain{k}.discrMethod = 'Analytic';
		% choose: 'Analytic', 'Numerical'
		%	Sets way of evaluation of integrals
		%	Analytic only for 'Leggett_hard' and 'LogZ'. Uses modified scheme by �itko
		%
		para.chain{k}.discretization = 'LogZ';
		% choose: 'LogZ','Linear'

		%%
		para.chain{k}.Lambda=1.2;						% Bath Discretization parameter
		para.chain{k}.z=1;                               % z-shift of bath; see Zitko 2009 - 10.1103/PhysRevB.79.085106
% 		para.chain{k}.L=0;                               % Length per bath; if L=0: choose optimal chain length according to para.precision;
% 		if L > 0								% chain length override
% 			para.chain{k}.L = L;
% 		end
		para.rescaling = 1;                     % rescale h1term, h2term for bosonchain with \lambda^{j-2}*h1term

	% Consistency checks
		if strcmp(para.chain{k}.discrMethod,'Analytic')
			assert(strcmp(para.chain{k}.spectralDensity,'Leggett_Hard'),'Analytic discretization only for Leggett with hard cutoff!');
		end
		if para.chain{k}.Lambda == 1
			para.chain{k}.discretization = 'Linear';
		elseif para.chain{k}.Lambda > 1
			para.chain{k}.discretization = 'LogZ';
		else
			error('VMPS:ModelDefinition','Lambda >= 1!')
		end
	else
		error('VMPS:ModelDefinition','You need to define para.chain.mapping!');
	end
end

%% Starting MPS Dimensions
D = 3;
dk = 20;
d_opt = 5;

if strcmp(para.model,'MLSpinBoson')     % definitions needed in SBM_genpara for spectral function & Wilson chain
    % Phonon Model Definition para.MLSB_mode:
    %   1:  from diagonalised, constant spacing Delta, predefined t=[t1 t2 t3 t4...]; energies symmetric about 0s
    %       Needs Define: MLSB_Ls, MLSB_Delta, MLSB_t,
    %       Automatically defined: SBM J(w)
    %   2:  Hamiltonian with rotational symmetry. Read in data from file.
    %       Needs Define: MLSBM_t, MLSB_system,
    %       Automatically defined: MLSB_Ls, Renger2012 J(w),
	if ~isempty(strfind(para.chain{1}.spectralDensity,'Leggett'))
		para.MLSB_mode = 1;
	elseif ~isempty(strfind(para.chain{1}.spectralDensity,'Renger'))
		para.MLSB_mode = 2;
	end
end

if isempty(strfind(para.model,'folded'))
	para.foldedChain = 0;                       % parameter to tell that Supersites for chain are used!
	para.M = 2;									% Number of 2-operator interaction terms per site in Hamiltonian.
else
	para.foldedChain = 1;
	para.M = 2*2;
	dk = dk^2;
end

para.spinposition=1;                            % This indicates all positions ~= bosonic! important for Vmat! The y chain is on the left and the z chain is on the right. (could be array !)
para.complex=0;                                 % set to 1 if any complex parameters are used.
para.resume=0;                                  % Read from saved results if available.
para.logging = 1;                               % Switch on logging and
parity = 0;										% 0: none; 1: odd, 2: even
para.precision = 5e-15;                         % was 5e-15; Determines chain length if L=0; Also E-convergence


%% %%%%%%% Calculate Wilson Chain parameters %%%%%%%%%%%%%%%%%%
% needed here: para.model, [para.MLSB_mode]
for k = 1:para.nEnvironments
	[para.chain{k}]=SBM_genpara(para.chain{k});               % only need alpha, s, Lambda. Returns epsilon and t of Wilson chain. Auto choose L if L == 0
	if (L == 0)
		L = para.chain{k}.L;                         % Take best L if not specially defined
		para.L = para.chain{k}.L;					 % not the best solution.
	end
end

%%
nc = para.nChains;											% = 1 for single chain models
para.D							= D*ones(1,L-1);			% Bond dimension; starting dimension is 2. para.D(L) is useless in this program
para.dk_start					= dk;						% local dimension per boson in bath. Will be increased effectively by oscillator shift.
para.dk							= para.dk_start*ones(nc,L);
para.dk(1,para.spinposition)	= 2;						% Impurity dimension
para.d_opt						= d_opt*ones(1,L);			% Dimension of first site is 2 (spin); Optimal Boson Basis dimension, was 16*ones
% TODO: d_opt = nc+1 x L if nc > 1; d_opt(i,:) addresses chain i, while
% d_opt(end,:) is the old OBB facing the MPS
para.d_opt(1,para.spinposition) = 2;						% Optimal Impurity dimension
para.eigs_tol					= 1e-8;
para.loopmax					= 50;
para.increasedk					= 0;						% Tells by how much dk should have been increased to achieve good sv in MPS. start with 0.

if strcmp(para.model,'SpinBoson2C')
	para.M = 4;
	para.dk(2,para.spinposition) = 1;		% non-existent singleton!
elseif strcmp(para.model,'SpinBoson3C')
	para.M = 6;
	para.dk(2,para.spinposition) = 1;		% non-existent singleton!
	para.dk(3,para.spinposition) = 1;		% non-existent singleton!
end

if strcmp(para.model,'2SpinPhononModel')
   %para.hz = para.hx;                  % use hz as level splitting. hx = hy = 0
   %para.hx = 0;
   para.hy = 0.1;                       % coupling between the 2 spin sites = J
   para.M = 4;                          % 4 terms in sum per site
   para.dk(1) = 4;                      % As kron(site1,site2), dim=4 on first site!
   para.d_opt(1) = 4;
   para.foldedChain=1;
end

%% Multi-Level Spin Boson Model for PPC

if strcmp(para.model,'MLSpinBoson')
%% Model Definition para.MLSB_mode:
    %   1:  from diagonalised, constant spacing Delta, predefined t=[t1 t2 t3 t4...]; energies symmetric about 0s
    %       Needs Define: MLSB_Ls, MLSB_Delta, MLSB_t,
    %       Automatically defined: SBM J(w)
    %   2:  Hamiltonian with rotational symmetry. Read in data from file.
    %       Needs Define: MLSBM_t, MLSB_system,
    %       Automatically defined: MLSB_Ls, Renger2012 J(w),
    %       Can be affected with static random disorder

% All modes:
% System-Bath coupling:
%    para.MLSB_t = [-1 -1 -1 1 1 1];                 % coupling of each level to bath from range [-1,1]. Used as para.t(1)*para.MLSB_t
    para.MLSB_p = period;                           % period of coupling for cosine
    para.MLSB_etaFactor = eta;                      % multiplier to increase system-bath coupling
    para.MLSB_t = @(x) para.MLSB_etaFactor.*exp(1i*2*pi/para.MLSB_p.*x);                % MLSB_t can be an anonymous function.
    para.MLSB_tOff = 1;                             % 0=diag coupling; 1 = 1. off-diagonal n,n+1 coupling

    if para.MLSB_mode == 1
        para.MLSB_Ls = 6;                           % number of levels in System, symmetric around E=0;
        para.MLSB_Delta = para.hx;                  % level spacing in already diagonalized model
    end

    if para.MLSB_mode == 2
        para.MLSB_system = 'RsMolischianumB850';
            % Choose: 'RsMolischianumB850', 'RsMolischianumB800B850'
        para.MLSB_Ls = length(Hamiltonian_PPC(para));
            % Inherently defined within model
    end

%% random static disorder:
    para.MLSB_staticDisorder = 1;
    para.MLSB_disSDV = 200/2.3548;         % standard deviation of disorder; FWHM = 200 cm^-1; FWHM = 2.3548*SDV;
    % Vector with absolute values of disorder in cm^-1:
    para.MLSB_disDiag = para.MLSB_disSDV.*randn(para.MLSB_Ls,1);        % set to 0 if no diag disorder wanted

%% All modes:
    para.dk(1) = para.MLSB_Ls;                  % also set dk properly
    para.d_opt(1) = para.MLSB_Ls;

end

if strcmp(para.model,'SpinBoson') || strcmp(para.model, 'SpinBoson2folded') || strcmp(para.model,'SpinBoson2C')|| strcmp(para.model,'SpinBoson3C')
%% Set-up parameters for specific ground state preparation!
    para.SpinBoson.GroundStateMode = 'artificial';
        % choose: 'decoupled', 'coupled', 'artificial';
		% -artificial does no optimization! this only sets up an artificial
		%		ground state with <n> = 0 on chain and InitialState 'sz'
    para.SpinBoson.InitialState = 'sz';
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
end

%%
para.SVDmethod	  = 'qr';					% 'qr': uses QR instead of SVD wherever possible; 'svd': use SVD always (slower) (Not working now)
% smallest SV shall lie between [max, min] otherwise truncate or expand
% the smaller the higher accuracy
para.svmaxtol	  = 1e-6;
para.svmintol	  = 1e-8;					% para.svmaxtol/2; %The lower limit for the smallest Vmat singular values.
para.adjust		  = 0;						% Initialize as 0, no need the edit. To adjust D. Is set = 1 in minimizeE.m
para.Dmin		  = 4;						% set a minimum Bond dimension.
para.dimlock	  = 0;						% set to 0 will change D and d_dop adaptively
para.minDimChange = 0.01;					% sets dimlock = 1 if relative dimension change < minDimChange. (larger makes less loops)

%% Parity settings
switch parity
    case 0
        para.parity='n';	% none
    case 1
        para.parity='o';	% odd
    case 2
        para.parity='e';	% even
end
para.spinbase='Z';
if para.parity~='n'
    para.Anzi=cell(para.L,1);
    para.Vnzi=cell(para.L,1);
    para.Anzilr=para.Anzi;
    para.Anzirl=para.Anzi;
end

%% %%%%%%%%%%%%%%%%%% OBB - Vmat related parameters %%%%%%%%%%%%%%%%%%%%
% Introduces Optimized Boson Basis (OBB)
para.useVmat=1;
if para.useVmat==0
    % then: d_opt = dk
    fprintf('Not using Vmat and OBB!\n');
    para.d_opt = para.dk;
     assert(para.dk_start==max(para.d_opt));
end
para.d_opt_min = 2;                                     % minimum d_opt dimension

%% %%%%%%%%%%%%%%%%%% dk Expansion - related parameters %%%%%%%
% only works together with OBB!
para.useDkExpand     = 0;		% Enable dk expansion, own algorithm. General switch

para.dkmax			 = 900;		% has to be square number for folded chains
para.expandprecision = 1e-5;	% unused?
para.hasexpanded     = 0;
% Method 1:
% needs para.useVmat = 1
para.useDkExpand1	 = 0;		% Expand dk if largest SV of site is smaller than this thershold. EMPIRICAL. Below this value, Vmat seems to need higher dk
para.expandBelowSV   = 0.995;
% Method 2:
% needs para.useVmat = 1
para.useDkExpand2	 = 1;       % Expand if wavefunction on site occupies the high energy dimensions
para.dkEx2_tail		 = 0.4;     % tail length [0 1] of occupation in Vmat to analyse
para.dkEx2_maxDev	 = 1.5;     % if std(log10(tail)) < maxDev --> no increase; Measures orders of magnitude in fluctuations of tail.
para.dkEx2_minExp	 = 13;      % if tail below this order than do not expand.
% Method 3:
% TODO:
%  analyse based on Amat (in case para.useVmat == 0)
%  Do the SVD in the local bond to create a Vmat and use it to estimate
%  tails.
if para.useDkExpand
	assert(para.useVmat == 1, 'dk expansion only works together with OBB Vmat');
end

%% %%%%%%%%%%%%%%%%%% Boson Shift - related parameters %%%%%%%%%%%%%%%%
% Introduces shift of bosonic oscillators to increase effective dk
para.useshift=0;
% only choose one of the following Methods
para.useFloShift = 0;                                   % shift all sites if trustsite > 0
para.useFloShift2 = 1;  para.FloShift2minSV = 0.995;    % shift everything if maximum Vmat SV fall below this value
para.useFloShift3 = 0;  para.FloShift3minMaxSV = 1;     % shift every 3rd loop. Shifts only sites where max SV worse than half of the worst SV
para.FloShift3loops = 5;                                % FloShift3 needs logging of results.Vmat_sv!
para.useChengShift=0;                                   % shifts sites < trustsite
para.useEveryShift=0;                                   % shift in every loop

para.shift=zeros(nc,para.L);							% separate shift for each chain!
para.relativeshift=zeros(nc,para.L);
para.relativeshiftprecision=0.01;        %When the relative shift is below this value then stop shifting
para=maxshift(para);

%% calculate shift analytically  (for independent sub-ohmic model)
%t = para.t; e = para.epsilon;
%A = gallery('tridiag',para.t(2:end),para.epsilon(1:end),para.t(2:end));    %creates tridiag for system A.(shiftVec) = (-t1*sigmaZ, 0,0,0,0,...)
%B = zeros(para.L-1,1);
%B(1) = -para.t(1)*1;
%para.shift =[0; A\B.*sqrt(2)]';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, name] = system('hostname');
para.hostname = strtrim(name);						% save hostname for later reference
para.version = 'v51';
if ~strcmp(computer,'PCWIN64')
	para.version = sprintf('%sTCM%s',para.version,para.hostname(3:end));
end
if strcmp(para.chain{1}.mapping,'OrthogonalPolynomials')
	para.version = ['OrthPol-',para.version];
elseif strcmp(para.chain{1}.discretization,'LogZ')
	para.version = ['LogZ-',para.version];
end

if para.chain{1}.s ~= 1
	para.version = sprintf('%s-s%g',para.version,para.s);
end

para.folder=sprintf([datestr(now,'yyyymmdd-HHMM'),'-%s-%s-alpha%.10gdelta%.10gepsilon%.10gdk%.10gD%.10gdopt%gL%d'],...
    para.model,para.version,alpha,delta,epsilon,dk,D,d_opt,L);
if strfind(para.model,'SpinBoson') && strcmp(para.SpinBoson.GroundStateMode,'artificial')
	para.folder = sprintf('%s-art-%s',para.folder,para.SpinBoson.InitialState);
end
para.filename=strcat(para.folder,'/results.mat');
if ~exist(para.filename,'file')
    mkdir(para.folder);
    para.savedexist=0;
else
    para.savedexist=1;
end

%% Start Calculation
[op,para]=genh1h2term(para);
[mps, Vmat,para,results,op] = minimizeE(op,para);

if ~isempty(strfind(para.model,'SpinBoson'))
%% Reset original parameters after specific ground state preparation!
    if strcmp(para.SpinBoson.GroundStateMode, 'decoupled')
		% restore coupling to chain
        para.t(1) = para.SpinBoson.t1;
    end
    if (strcmp(para.SpinBoson.InitialState, 'sx') || strcmp(para.SpinBoson.InitialState, 'sz')) && ~strcmp(para.SpinBoson.GroundStateMode, 'artificial')
        para.hx = para.SpinBoson.hx;
        para.hz = para.SpinBoson.hz;
    end
    [op,para]   = genh1h2term(para);                % restore operators!
    [op]        = initstorage(mps, Vmat, op,para);

end


save(para.filename,'para','Vmat','mps','results','op','-v7.3');

%% Calculate some Results
results.nx         = getObservable({'occupation'},mps,Vmat,para);
% results.bosonshift = getObservable({'shift'},mps,Vmat,para);

if strcmp(para.model,'SpinBoson') || strcmp(para.model,'SpinBoson2folded')
    results.spin   = getObservable({'spin'},mps,Vmat,para);
end

if strcmp(para.model,'MLSpinBoson')
    results.participation = getObservable({'participation'},mps,Vmat,para);
    results.tunnelEnergy  = getObservable({'tunnelenergy',op},mps,Vmat,para);
end
%%
results.time = toc(starttime)
save(para.filename,'para','Vmat','mps','results','op','-v7.3');

fileName = para.filename;
end
