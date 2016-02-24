function tdvp_1site_tree(treeMPS,para,results,tresults)
%% VMPS Time-dependent variational principle
%   implementation following Haegeman et al. 2014, arXiv:1408.5056
%   using Expokit for expv()
%
%   Expands OBB on sitej upon entering mps{sitej}, before SVD to Vmat
%   Truncates OBB on sitej after focusing the next mps site
%
%
%   To continue a completed calculation use:
%       inarg: tmps instead of mps, tVmat instead of Vmat. reorganization
%              will be done in code under 0.
%
%	To resume an incomplete calculation use:
%		inarg: tmps(end)-> mps, tVmat(end) -> Vmat
%
%   Created by Florian Schroeder 07/10/2014
%
%	Modified 22/02/2016
%		- Restructured Variables to support TreeMPS
%

if nargin == 0
	createTestRun();
end

%% 0. Parameter settings:
if ~isfield(para,'tdvp')
    error('VMPS:tdvp_1site_tree:MissingParameters','Please define all necessary tdvp variables in para!');
end

% open file for writing results
% if storeMPS then para.tdvp.filename will contain tMPS
assert(isfield(para.tdvp,'filename'),'Need definition of para.tdvp.filename to save results!');
assert(isfield(para.tdvp,'filenameSmall'),'Need definition of para.tdvp.filenameSmall to save results!');

if isfield(para.tdvp,'storeMPS') && para.tdvp.storeMPS == 1
	outFile = matfile(para.tdvp.filenameMPS,'Writable',true);
	% this does not yet create file.
	% File will be created upon first writing!
end

resumeTDVP = 0;					% 0 = start new, 1 = resume aborted, 2 = continue completed

if ~isempty(tresults) && para.tdvp.resume == 1 && exist(para.tdvp.filename,'file')
	%% must be resumed since no writing has happend since start!
	% detect startTime:
	%  - outFile exists, and is already synced with tresults from resume
	%  - only have tresults
% 	try
% 		load(para.tdvp.filename);	% all the rest
% 	catch
% 		disp('MAT-file corrupt, load BAK');
% 		load([para.tdvp.filename(1:end-4),'.bak'],'-mat');		% use backup if file corrupted
% 	end
	if isfield(tresults,'lastIdx')
		% there must be preexisting calculation -> resume!
		startSize = para.timeslice+1;
		resumeTDVP = 1;
	else
		startSize = 1;				% fallback value. Perhaps change?
		resumeTDVP = 0;
		clear('tresults');			% throw away old results -> Overwrite!
	end
elseif ~isempty(tresults)
	% NEEDS WORK!! TODO
	% detect continuation
	if isempty(strfind(para.tdvp.fromFilename,'results.mat'))
		resumeTDVP = 2;
		%load(para.tdvp.fromFilename, 'tresults');							% if error, load .bak file!
		% tresults will be expanded in calTimeObservables!

% 		if ~exist('tresults','var')
% 			error('tresults were not in fromFile!');
% 		end
		startSize = para.timeslice+1;
		fprintf('\nSet resumeTDVP = 2, startSize = %d', startSize);
	end
end

%% 1. sweep: condition the ground state MPS and prepare hamiltonian operator terms H(n) and K(n)
    % is already done in op.opstorage and op.hlrstorage
    %
    % also prepare other stuff
if resumeTDVP == 0
	if exist('outFile','var')
		% initialize tMPS and calc observables at the end
		fprintf('Initialize tMPS\n');
		error('VMPS:tdvp_1site_tree:NotImplemented','This needs implementation.')
% 		outFile.tmps(1,para.L) = {[]};
% 		outFile.tVmat(1,para.L) = {[]};
% 		outFile.tmps(1,:) = mps;
% 		outFile.tVmat(1,:) = Vmat;
	end
	% calc tresults now to initialize tresults
% 	tresults = calTimeObservables(treeMPS,para);				% TODO: this needs overloading..
	if para.logging
% 		results = initresultsTDVP(para, results);				% TODO: this needs Work
	end
	para.trustsite(end) = para.L;       % necessary to enable several truncation / expansion mechanisms
	% For fresh starts, define number of timeslices
	para.tdvp.slices = 1:(length(para.tdvp.t)-1);
	para.tdvp.calcTime = zeros(length(para.tdvp.t)-1,1);		% logs the elapsed time after each slice in h
elseif resumeTDVP == 1
	fprintf('Resuming aborted TDVP\n')
	para.tdvp.slices = startSize:(length(para.tdvp.t)-1);
	if ~isfield(results,'tdvp')
		% results was not stored properly before -> reinitialize
		if para.logging
% 			results = initresultsTDVP(para, results);			% TODO!!!
			if para.tdvp.expandOBB
% 				results.tdvp.d_opt(startSize,:)   = para.d_opt;		% TODO!!
			end
			if para.tdvp.truncateExpandBonds
% 				results.tdvp.D(startSize,:)       = para.D;		% TODO!!
			end
		end
	end
elseif resumeTDVP == 2
	%% If continue previously completed TDVP
	fprintf('Continuing previous TDVP\n');
	para.tdvp.slices = startSize:(length(para.tdvp.t)-1);
	if exist('outFile','var');								% TODO!!
% 		outFile.tmps = mps;
% 		outFile.tVmat = Vmat;
% 		outFile.tmps(length(para.tdvp.t),para.L) = {[]};
% 		outFile.tVmat(length(para.tdvp.t),para.L) = {[]};
	end
	% increase size of results objects:
	if para.logging
		n = length(para.tdvp.t);
		if para.tdvp.expandOBB
% 			results.tdvp.d_opt(n,end)   = 0;						% TODO!!
		end
		if para.tdvp.truncateExpandBonds
% 			results.tdvp.D(n,end)       = 0;						% TODO!!
		end
% 		results.tdvp.expError(n-1,1)    = 0;						% TODO!!
% 		results.tdvp.EvaluesLog(n-1,1)  = 0;						% TODO!!
		if para.tdvp.logSV
% 			results.tdvp.Vmat_sv{n,end} = cell(n,para.L);						% TODO!!
% 			results.tdvp.Amat_sv{n,end} = cell(n,para.L-1);						% TODO!!
		else
% 			results.tdvp.Vmat_vNE(n,end)= 0;						% TODO!!
% 			results.tdvp.Amat_vNE(n,end)= 0;						% TODO!!
		end
	end
end

fprintf('Saving.');
save(para.tdvp.filenameSmall,'para','tresults','-v7.3');
if para.tdvp.serialize
	hlp_save(para.tdvp.filename,para,results,tresults,treeMPS)
else
	save(para.tdvp.filename,'para','results','tresults','treeMPS','-v7.3');
	fprintf('.\n');
end
%% Asserts to do:

%% 2. time sweep

for timeslice = para.tdvp.slices
    para.sweepto = 'r'; treeMPS.sweepto = 'r';
	para.timeslice = timeslice;
% 	para.tdvp.deltaT = para.tdvp.t(timeslice+1)-para.tdvp.t(timeslice);
    fprintf('t = %g\n', para.tdvp.t(timeslice+1));
    % sweep l->r and time evolve each site
	
	%% check whether treeMPS is Tree, Star or Chain
	%	Pro-Recursive: No need for call stack, automatic assignment after calculation
	%	Con-Recursive: Function call overhead. copy on write
	%	Pro-Iterative: No Function call overhead
	%	Con-Iterative: Call Stack needed. Saving to parent problematic, as not doable by reference. 
	
	% Do recursion now!
	if treeMPS.height > 0
		% treeMPS is Tree (or Star)
		treeMPS = tdvp_1site_evolveTree(treeMPS,para);
	else
		% treeMPS is a Leaf/Chain
		treeMPS = tdvp_1site_evolveChain(treeMPS,para);
	end
	
	%% Perform Site 1 TDVP
	% only use expvCustom due to complicated contractions!
	fprintf('0');
	para.sitej = 1;
	mps = tdvp_1site_evolveSystem(mps,Vmat,para,op,results);			% nested function for now!
	fprintf('-');
% 	Cn = 0;	% declare for scope


	%% run a HOSVD to focus each chain
	for mc = 1:para.nChains
	%% Now normal sweep through chains. Redefine op to bring onto single-chain form!
		fprintf('\nChain %d:\n',mc);
		para.currentChain = mc; para.sitej = 1;
		% SVD to get A_C(mc)
		[mps{1}, Cn] = prepare_Tens(mps{1}, para);		% mps{1} can be seen as Vtens core

		% create site 1 OBB terms h1j/h2jOBB
		[op] = H_Eff([]    ,Vmat{1}, 'A'    , op, para);
		% calculate op.Hleft, op.Opleft for the next chain!
		[op] = H_Eff(mps{1}, []    , 'ST-CA', op, para);
			%also need H/Opright!
		% prepare to form single chain
		[mpsChain, VmatChain, paraChain, opChain, resultsChain] = prepare_ChainOp(mps, Vmat, para, op, results);

		% evolve AC back in time
		paraChain.tdvp.expvCustomNow = 1;		% avoid need for Hn
		paraChain.sweepto = 'r';

		% declare site as 0 to achieve contraction into mpsChain{1}
		[mpsChain, ~, paraChain, resultsChain] = tdvp_1site_evolveKn(mpsChain,[],paraChain,resultsChain,op,0,Cn,[]);
		% now focused on mps{2}{mc} == mpsChain{1}
		% start the standard chain-sweep!

		%%
		
	
		% save first special Bond info by Hand
		results.Amat_sv{mc,1}  = sv;
		results.Amat_vNE(mc,1) = vNE;		% from last prepare_onesite_truncate

		%% finished with l<-r sweep, focused on CA(1) after time-evolution of chain mc
		para.D(mc,1) = size(Cn,2);
		% Hlrstorage{1} already updated loop in updateop!
		op.chain(mc).Hlrstorage = opChain.Hlrstorage;
		op.chain(mc).Opstorage  = opChain.Opstorage;
		% just for clarity:
		op.h1j = []; op.h2j = [];
		% Contract Cn into mps{1} and update all other StarMPS Variables
		NC = para.nChains;
		mps{1} = contracttensors(mps{1}, NC+2, mc+1, Cn, 2, 1);
		mps{1} = permute(mps{1}, [1:mc, NC+2, mc+1:NC+1]);
		[mps,Vmat,para,results] = copyChainToStar(mps,mpsChain,Vmat,VmatChain,para,paraChain,results,resultsChain);
	end
	% Evolve mps{1} by dt/2 to finish one timesweep!
	para.sitej = 1;
	mps = tdvp_1site_evolveSystem(mps,Vmat,para,op,results);		% nested function for now!

	% normalise the state
	d = size(mps{1}); para.sweepto = 'l';
	[mps{1}, Cn] = prepare_onesite(reshape(mps{1},1,[],d(end)),para,1);
	mps{1} = reshape(mps{1},d);


	fprintf('\n');
	if paraChain.tdvp.truncateExpandBonds
		fprintf('para.D:\n');
		out = strrep(mat2str(para.D),';','\n');
		fprintf([out(2:end-1),'\n']);
		out = strrep(mat2str(cellfun(@(x) x(find(x,1,'last')),results.Amat_sv),2),';','\n');
		fprintf([out(2:end-1),'\n']);
	end
	if paraChain.tdvp.expandOBB
		fprintf('para.d_opt:\n');
		out = strrep(mat2str(para.d_opt),';','\n');
		fprintf([out(2:end-1),'\n']);
		out = strrep(mat2str(cellfun(@(x) -log10(x(end)),reshape(results.Vmat_sv(~cellfun('isempty',results.Vmat_sv)),para.nChains,[])),2),';','\n');
		fprintf([out(2:end-1),'\n']);
	end


	% report time estimates
    fprintf('Norm: %g\n', Cn);
	if sitej == 1
		results.tdvp.EvaluesLog(timeslice) = getObservable({'energy',op},mps,Vmat,para);
	end
	completePercent = round(timeslice./(length(para.tdvp.t)-1)*1000)./10;
	para.tdvp.calcTime(timeslice) = toc(para.tdvp.starttime)./3600;
	n = 1:timeslice;
	n = n(para.tdvp.calcTime(n,1)>0);		% compensates for zeros due to resume of aborted TDVP
	hoursTotal = fnval(fnxtr(csaps([0,n],[0;para.tdvp.calcTime(n,1)])), length(para.tdvp.t)-1);
	if log10(abs(hoursTotal-para.tdvp.calcTime(timeslice))) <= 0
		fprintf('Completed: %5.3g%%, t = %4.3g, Time elapsed: %4.2gh, Time left: %4.2gm, Time total: %4.2gh\n',...
			completePercent, para.tdvp.t(timeslice+1), para.tdvp.calcTime(timeslice),...
			(hoursTotal-para.tdvp.calcTime(timeslice))*60, hoursTotal);
	else
		fprintf('Completed: %5.3g%%, t = %4.3g, Time elapsed: %4.2gh, Time left: %4.2gh, Time total: %4.2gh\n',...
			completePercent, para.tdvp.t(timeslice+1), para.tdvp.calcTime(timeslice),...
			hoursTotal-para.tdvp.calcTime(timeslice), hoursTotal);
	end
    %% save tmps and tVmat and log parameters
	if exist('outFile','var')
		fprintf('Saving tMPS\n');
		outFile.tmps(timeslice+1,:)  = mps;					% writes to File
		outFile.tVmat(timeslice+1,:) = Vmat;				% writes to File
		outFile.currentSize = timeslice+1;
	end
	if para.logging && mod(timeslice, round(para.tdvp.extractObsInterval/para.tdvp.deltaT))==0			% at extractObsInterval
		n = 1+tresults.lastIdx;
		if para.tdvp.logSV
			if size(results.Vmat_sv,1) == para.L
				results.tdvp.Vmat_sv(n,:,:)	 = results.Vmat_sv;
			else
				results.tdvp.Vmat_sv(n,:,:)	 = results.Vmat_sv.';
			end
			results.tdvp.Amat_sv(n,:) = results.Amat_sv;
		else
			results.tdvp.Vmat_vNE(n,:,:)= results.Vmat_vNE.';
			results.tdvp.Amat_vNE(n,:,:)= results.Amat_vNE.';
		end
		if para.tdvp.expandOBB
	        results.tdvp.d_opt(n,:)   = reshape(para.d_opt.',1,[]) - sum(results.tdvp.d_opt);	% sum in 1st Dim
		end
		if para.tdvp.truncateExpandBonds
	        results.tdvp.D(n,:)       = reshape(para.D.',1,[]) - sum(results.tdvp.D);
		end
		if para.tdvp.useDkExpand
			results.tdvp.dk(n,:)      = reshape(para.dk.',1,[]) - sum(results.tdvp.dk);
		end
	end
	if ~exist('tresults','var')
		% calculate everything up to now!
		if exist('outFile','var')
			tresults = calTimeObservables(outFile.tmps(1:outFile.currentSize,:),outFile.tVmat(1:outFile.currentSize,:),para);
		else
			% Should not happen!
			error('How did you get here??');
			%tresults = calTimeObservables(mps,Vmat,para);
		end
	elseif isfield(tresults,'lastIdx') % && tresults.lastIdx == timeslice not applicable anymore for extractObsInterval
		% only calculate the current slice
		tresults = calTimeObservables(mps,Vmat,para,tresults);
	end


	%% Additional saving
	if mod(timeslice,para.tdvp.saveInterval) == 0 || timeslice ==  para.tdvp.slices(end)
		%% backup results, op and para less often than tmps and tVmat
		results.tdvp.time = toc(para.tdvp.starttime);
		fprintf('Saving results.. ')
		copyfile(para.tdvp.filename,[para.tdvp.filename(1:end-4),'.bak']);			% in case saving gets interrupted
		save(para.tdvp.filenameSmall,'para','tresults','-v7.3');
		if para.tdvp.serialize
			hlp_save(para.tdvp.filename,para,op,results,tresults,mps,Vmat)
		else
			save(para.tdvp.filename,'para','op','results','tresults','mps','Vmat','-v7.3');
			fprintf('.\n');
		end
	end
end

delete([para.tdvp.filename(1:end-4),'.bak']);			% get rid of bak file

    function truncateOBB(sitej)
        %% Truncate OBB of mps{sitej}
        % code taken from adjustdopt.m
        % only apply if also expanded OBB before
		if para.nChains > 1, return, end
        discarddims     = find(results.Vmat_sv{sitej} < para.svmintol);
        if ~isempty(discarddims) && para.tdvp.expandOBB
            if results.Vmat_sv{sitej}(discarddims(1)-1) > para.svmaxtol     % if next highest not-discarded element too large (causing expansion in next sweep)
                discarddims = discarddims(2:end);                           % always keep one element < svmaxtol
            end
            % if Dif > 0: would remove too many dims; if Dif < 0, no
            % problem so set Dif = 0
            difference = para.d_opt_min + length(discarddims) - para.d_opt(sitej);
            if difference <= 0      % discarddims does not violate d_opt_min
                difference = 0;
            end                     % else: discarddims would remove too much by amount = diff;
%                 dispif('remove dims in d_opt',para.logging)
            Vmat{sitej}(:,discarddims(difference+1:end))=[];                          % only cut the end
            mps{sitej}(:,:,discarddims(difference+1:end))=[];
            para.d_opt(sitej)=size(Vmat{sitej},2);
            results.Vmat_sv{sitej}(discarddims(difference+1:end))=[];
        end
	end

	function [Vmat,para,paraC,results,op,opC] = truncateExpandDk(Vmat,para,paraC,results,op,opC,s)
			%% function [Vmat, para, op] = truncateExpandDk(Vmat,para,op)
			%	truncates or expands the local dimension dk based on the OBB usage
			%	the vectors in Vmat and their SV are the indicators
			%	for now only for standard single-chain VMPS
			%	Input: Vmat = Vmat{sitej}, single-site Vmat
			if para.useVtens || ~para.useDkExpand || ~para.useVmat
				return;			% Do nothing
			end
			% s = sitej;
			nc = para.currentChain;
			M = para.M/para.nChains;

			%% Determine whether to truncate or to expand: Copied from adjustdopt.m
			adddim = estimateDkExpand(Vmat,results.Vmat_sv{1,s}, para);
			if adddim > 0
			%% Apply the expansion
				[dk,dOBB]   = size(Vmat);
				paraC.dk(s)  = dk+adddim;						% operators will be expanded in genh1h2term?
				para.dk(nc,s+1) = paraC.dk(s);				% needs to be updated now!
				addmat  = zeros(adddim,dOBB);
	% 			Vmat = cat(1,addmat,Vmat);			% N:1 ordering
				Vmat = cat(1,Vmat,addmat);			% 1:N ordering
				paraC.hasexpanded = 1;
				paraC.increasedk  = 0;							% reset this value
	% 			dispif('Increased dk',para.logging)
			end
			% expand the operators of the Full StarMPS:
			op = genh1h2term_onesite(para,op,s+1);			% additional s+1, since in StarMPS picture
			% copy out the chain operators
			%   although they might not be needed anymore
			opC.h1term(s)     = op.h1term(nc, s+1);
			opC.h2term(:,:,s) = op.h2term( M*(nc-1)+(1:M), :, s+1, nc);

			% since only zeros were added in Vmat, no further update of h1j/OBB etc is needed.
	end

	function results = initresultsTDVP(para, results)
		% save Dimension Log only as difference! reconstruct with cumsum(A)
% 		fprintf('Initialize results.tdvp\n');
		n = round(para.tdvp.tmax/para.tdvp.extractObsInterval) +1;
		NC = para.nChains; L = para.L;
		if para.tdvp.logSV
			if para.useVtens
				results.tdvp.Vmat_sv{n,L,NC+1}	= [];
			elseif para.useStarMPS
				results.tdvp.Vmat_sv{n,L,NC}	= [];
			end
			if size(results.Vmat_sv,1) == L
				results.tdvp.Vmat_sv(1,:,:)	 = results.Vmat_sv;
			else
				results.tdvp.Vmat_sv(1,:,:)	 = results.Vmat_sv.';
			end
			results.tdvp.Amat_sv{n,L-1} = [];
			results.tdvp.Amat_sv(1,:)		 = results.Amat_sv;
		else
			if para.useVtens
				results.tdvp.Vmat_vNE(n,L, NC+1) = 0;
				results.tdvp.Amat_vNE(n,L-1)	 = 0;
			elseif para.useStarMPS
				results.tdvp.Vmat_vNE(n,L, NC)   = 0;
				if L ~= 1
					results.tdvp.Amat_vNE(n,L-1, NC) = 0;
				end
			end
			results.tdvp.Vmat_vNE(1,:,:)	= results.Vmat_vNE.';
			results.tdvp.Amat_vNE(1,:,:)	= results.Amat_vNE.';
		end
		if para.tdvp.expandOBB
			if para.useVtens
				results.tdvp.d_opt              = sparse(n,L*(NC+1));
			else
				results.tdvp.d_opt              = sparse(n,L*NC);
			end
			results.tdvp.d_opt(1,:)			= reshape(para.d_opt.',1,[]);
		end
		if para.tdvp.truncateExpandBonds
			if para.useStarMPS
				results.tdvp.D				= sparse(n,(L-1)*NC);
			else
				results.tdvp.D				= sparse(n,L-1);
			end
			results.tdvp.D(1,:)				= reshape(para.D.',1,[]);
		end
		if para.tdvp.useDkExpand
			if para.useStarMPS
				results.tdvp.dk				= sparse(n,L*NC);
			else
				results.tdvp.dk				= sparse(n,L);
			end
			results.tdvp.dk(1,:)			= reshape(para.dk.',1,[]);
		end
%		results.tdvp.dk 		  = sparse(n,para.L);
%		results.tdvp.dk(1,:)      = para.dk;
		results.tdvp.expvTime     = [];
% 		results.tdvp.expError	  = zeros(length(para.tdvp.t)-1, (para.L-1)*4*2+3);
		results.tdvp.expError	  = zeros(length(para.tdvp.t)-1, 1);		% otherwise file gets too large
		results.tdvp.EvaluesLog	  = zeros(length(para.tdvp.t)-1, 1);
	end

	function [mpsChain, VmatChain, paraChain, opChain, rChain] = prepare_ChainOp(mps, Vmat, para, op, results)
		%% prepares mps, Vmat and op to use them in single-chain code
		% focus should be already on mps{2}{currentChain}
		% Created 08/10/2015 by FS

		nc = para.currentChain;
		L = para.chain{nc}.L;
		M = para.M/para.nChains;

		mpsChain    = cellfun(@(x) x{nc},mps(2:L),'UniformOutput',false);
		VmatChain   = cellfun(@(x) x{nc},Vmat(2:L),'UniformOutput',false);

		paraChain   = para;
		paraChain.L = L-1;
		paraChain.M = M;
		paraChain.nChains      = 1;
		paraChain.spinposition = [];					% there is no spin on subchain!
		paraChain.dk           = para.dk(nc,2:L);
		paraChain.d_opt        = para.d_opt(nc,2:L);
		paraChain.D            = para.D(nc,2:L-1);

		if ~isempty(results)
			rChain = initresults(paraChain);				% without any previous settings
			rChain = initresultsTDVP(paraChain, rChain);
		end
		% Copy operators
		opChain.h1term = cell(1,L-1);
		opChain.h2term = cell(M, 2,L-1);

		opChain.h1term = op.h1term(nc, 2:L);
		opChain.h2term = op.h2term( M*(nc-1)+(1:M), :, 2:L, nc);

		opChain = initstorage(mpsChain,VmatChain,opChain,paraChain);

		% need to find these terms: Take from H_Eff ??
		opChain.Hlrstorage{1, 1} = op.chain(nc).Hleft;				% H left non-interacting, from: H_Eff 'ST-CA'
		for m = 1:M
			opChain.Opstorage{m,1,1} = op.chain(nc).Opleft{m};		% H left interacting
		end

	end

	function createTestRun()
		%% function createTestRun()
		%	
		%	creates a testcase for debugging and development. Will be called if no arguments are given.
		%	
		Lchains = 5;
		D = 10;			% one BondDim for all edges
		d_opt = 10;
		dk = 20;
		
		% Define tree structure in para:
		para.logging = 1;
		para.L = 6;								% tree nodes count all into site 1
		para.useTreeMPS = 1;
		para.useStarMPS = 0;
		para.useVmat = 1;
		para.trustsite = para.L;
		para.model = 'testTree';
		para.nChains = 4;						% = external nodes = leaves
		para.nNodes = 2;						% = internal nodes = nodes
		para.alpha = 0.1;
		para.hx = 0.1;
		para.hz = 0;
		para.M = 2;
		para.spinbase = 'Z';
		para.parity = 'n';
% 		para.rescaling = 0;
		para.treeMPS.height = 2;
		para.treeMPS.maxDegree = [3,2]; % size(treeMPS.L);		% for each level
		para.treeMPS.leafIdx = num2cell(ones(para.nChains,para.treeMPS.height));				% indices of leaves in para; take -1 to get treeIdx
		para.treeMPS.leafIdx{1,1} = 1+1;														% map from chain number to leaf index in para
		para.treeMPS.leafIdx{2,1} = 2+1;
		para.treeMPS.leafIdx(3,1:2) = num2cell([3,1]+1);
		para.treeMPS.leafIdx(4,1:2) = num2cell([3,2]+1);
		para.treeMPS.nodeIdx = num2cell(ones(para.nNodes,para.treeMPS.height));					% mpas node number -> nodeIdx
		para.treeMPS.nodeIdx{1,1} = 0+1;
		para.treeMPS.nodeIdx{2,1} = 3+1;
		para.treeMPS.chainIdx = {};																% maps from leafIdx -> chain number
		for ii = 1:size(para.treeMPS.leafIdx,1)
			para.treeMPS.chainIdx{para.treeMPS.leafIdx{ii,:}} = ii;
		end
		para.chain{1}.alpha   = 0.1;
		para.chain{1}.L       = Lchains;
		para.chain{1}.s       = 1;
		para.chain{1}.epsilon = ones(Lchains,1)*0.5;
		para.chain{1}.t       = ones(Lchains,1)*0.25;
		for ii = 1:para.nChains
			para.chain{ii} = para.chain{1};
		end
		
		% Define each struct of chainMPS and initialise random MPS
		structChainMPS(para.nChains) = struct();		% smallest subunit
		for ii = 1:para.nChains
			Lc = para.chain{ii}.L;
			chainMPS = cell(1,Lc);
			chainVmat = cell(1,Lc);
			for jj = 1:Lc
				if jj == Lc
					chainMPS{jj} = randn(D,1,d_opt);			% end-of-chain -> Dr = 1
				else
					chainMPS{jj} = randn(D,D,d_opt);
				end
				chainVmat{jj} = sparse(1:d_opt,1:d_opt,1,dk,d_opt);
			end
			% Build Struct Chain MPS
			structChainMPS(ii).mps = chainMPS;
			structChainMPS(ii).Vmat = chainVmat;
			structChainMPS(ii).height = 0;						% this struct is a leaf
			structChainMPS(ii).degree = 0;						% leaf
			structChainMPS(ii).child = [];						% contains all children nodes
			structChainMPS(ii).currentChild = [];				% leaf: [], node: currentChain;
			structChainMPS(ii).L = getTreeLength(structChainMPS(ii));		% StarMPS can have array L
			structChainMPS(ii).D = [ones(1,Lc)*D,1];		% first is Dl of chain, last is Dr of chain, here horizontal!
			structChainMPS(ii).d_opt = ones(1,Lc)*d_opt;
			structChainMPS(ii).dk = ones(1,Lc)*dk;
			structChainMPS(ii).M = para.M;
			structChainMPS(ii).shift = zeros(1,Lc);
			structChainMPS(ii).useVmat = 1;
			structChainMPS(ii).useVtens = 0;
			structChainMPS(ii).useStarMPS = 0;
			structChainMPS(ii).useTreeMPS = 0;
			structChainMPS(ii).nChains = 1;
			structChainMPS(ii).sweepto = 'l';					% 'l' for initial preparation
			structChainMPS(ii).sitej = 1;						% currently focused site
			structChainMPS(ii).isRoot = 0;
			structChainMPS(ii).treeIdx = [para.treeMPS.leafIdx{ii,:}]-1;
			structChainMPS(ii).BondCenter = [];
			pIdx = para.treeMPS.leafIdx(ii,:);
			para.dk{pIdx{:}} = structChainMPS(ii).dk;
			para.D{pIdx{:}} = structChainMPS(ii).D;
			para.d_opt{pIdx{:}} = structChainMPS(ii).d_opt;
			para.shift{pIdx{:}} = structChainMPS(ii).shift;
			para.treeMPS.L(pIdx{:}) = structChainMPS(ii).L;
			
			% generate Hamiltonian terms:
			structChainMPS(ii).op.h1term = cell(1,Lc);
			structChainMPS(ii).op.h2term = cell(para.M,2,Lc);
			for jj = 1:Lc
				temp = genh1h2term_onesite(para,structChainMPS(ii).treeIdx,jj);
				structChainMPS(ii).op.h1term(jj)     = temp.h1term;
				structChainMPS(ii).op.h2term(:,:,jj) = temp.h2term;
			end
		end
		
		% Construct tree node 2
		ii = 2;
		pIdx = para.treeMPS.nodeIdx(ii,:);
		dk = 2; d_opt = 2;
		structStarMPS = structChainMPS(1);					% make a copy to inherit fields. Only modify difference
		structStarMPS.treeIdx = [para.treeMPS.nodeIdx{ii,:}]-1;
% 		[para.treeMPS.nodeIdx{ii,structStarMPS.treeIdx ~= 0}]
		structStarMPS.mps = {randn(D,D,D,d_opt)};		% Dl, Dc1, Dc2, dOBB
		structStarMPS.Vmat = {sparse(1:d_opt,1:d_opt,1,dk,d_opt)};
		structStarMPS.child = [structChainMPS(3);structChainMPS(4)];
		structStarMPS.currentChild = 1;
		structStarMPS.height = max([structStarMPS.child.height])+1;
		structStarMPS.degree = 2;
		structStarMPS.nChains = 2;						% equals degree for nodes
		structStarMPS.L = [];
		structStarMPS.L = getTreeLength(structStarMPS);		% StarMPS can have array L
		structStarMPS.d_opt = d_opt;
		structStarMPS.dk = dk;
		structStarMPS.D = ones(3,1)*D;					% [Dl;Dc1;Dc2]  here vertical!
		structStarMPS.shift = 0;
		structStarMPS.useStarMPS = 1;					% star is subset of tree!
		structStarMPS.useTreeMPS = 1;
		para.dk{pIdx{:}} = structStarMPS.dk;
		para.D{pIdx{:}} = structStarMPS.D;
		para.d_opt{pIdx{:}} = structStarMPS.d_opt;
		para.shift{pIdx{:}} = structStarMPS.shift;
		para.treeMPS.L(pIdx{:}) = 1;
		% generate Hamiltonian terms:
		structStarMPS.op = genh1h2term_onesite(para,structStarMPS.treeIdx,1);

		ii = 1;
		pIdx = para.treeMPS.nodeIdx(ii,:);
		dk = 2; d_opt = 2;
		treeMPS = structStarMPS;							% make a copy to inherit fields. Only modify difference
		treeMPS.mps = {randn(1,D,D,D,d_opt)};
		treeMPS.Vmat = {sparse(1:d_opt,1:d_opt,1,dk,d_opt)};
		treeMPS.child = [structChainMPS(1); structChainMPS(2); structStarMPS];
		treeMPS.currentChild = 1;
		treeMPS.height = max([treeMPS.child.height])+1;
		treeMPS.degree = 3;
		treeMPS.nChains = 2;						% equals degree for nodes
		treeMPS.L = [];
		treeMPS.L = getTreeLength(treeMPS);				% TreeMPS needs cell L for now..
		treeMPS.d_opt = d_opt;
		treeMPS.dk = dk;
		treeMPS.D = [1;ones(3,1)*D];					% [Dl;Dc1;Dc2;Dc3]  here vertical! Leading singleton -> root of tree!
		treeMPS.useStarMPS = 0;					% star is subset of tree!
		treeMPS.useTreeMPS = 1;
		treeMPS.isRoot = 1;
		treeMPS.treeIdx = [para.treeMPS.nodeIdx{ii,:}]-1;
		treeMPS.op = genh1h2term_onesite(para,treeMPS.treeIdx,1);
		para.treeMPS.L(pIdx{:}) = 1;
		
% 		para.treeMPS.L = treeMPS.L;		% not necessary anymore!
		
		para.tdvp.filename = 'test.mat';
		para.tdvp.filenameSmall = 'test-small.mat';
		para.tdvp.StoreMPS = 0;
		para.tdvp.resume = 0;
		para.tdvp.tmax = 3;
		para.tdvp.deltaT = 1;
		para.tdvp.t = 0:para.tdvp.deltaT:para.tdvp.tmax;
		para.tdvp.serialize = 0;
		para.tdvp.HEffSplitIsometry = 0;
		
		tresults = [];
		results = [];
		
		%% Do one prepare sweep to bring MPS onto right-canonical form
		[treeMPS,~,para] = prepare(treeMPS,[],para);
		
		%% Calculate effective Hamiltonians
		treeMPS = initstorage(treeMPS,[],[],para);		% new call pattern for treeMPS
		
	end
end

function [mps,Vmat,para,results] = copyChainToStar(mps,mpsC,Vmat,VmatC,para,paraC,results,resultsC)
% only use this after sweep along a single chain
nc = para.currentChain;
NC = para.nChains;
L  = para.chain{nc}.L;

% Copy the MPS from Chain to Star!
for ii = 1:length(mpsC)
	mps{ii+1}{nc} = mpsC{ii};
	Vmat{ii+1}{nc} = VmatC{ii};
end

para.dk(nc,2:L)    = paraC.dk;
para.d_opt(nc,2:L) = paraC.d_opt;
para.D(nc,2:L-1)   = paraC.D;

results.Vmat_sv(nc,2:L)    = resultsC.Vmat_sv;
results.Amat_sv(nc,2:L-1)  = resultsC.Amat_sv;
results.Vmat_vNE(nc,2:L)   = resultsC.Vmat_vNE;
results.Amat_vNE(nc,2:L-1) = resultsC.Amat_vNE;

end

function [Atens, CA] = prepare_Tens(Atens, para)
	%% Input: focused AS = mps{1} for star-MPS
	%  sets focus onto Vtens{mc}
	%  no truncation for now!
	%
	mc = para.currentChain;								% the chain to focus on
	d  = size(Atens);

	[AS]     = tensShape(Atens, 'unfold', mc+1, d);		% D(mc) x rest

	[AS, CA] = prepare_onesiteVmat(AS.', para);			% need for transpose since m x n, m > n input needed

	Atens	 = tensShape(AS.', 'fold', mc+1, d);

end



function treeMPS = tdvp_1site_evolveTree(treeMPS,para)
%% function treeMPS = tdvp_1site_evolveTree(treeMPS,para)
%
%	evolves node and calls evolution of children recursively
%	treeMPS is not shared with workspace of tdvp_1site_tree. Only shared with tdvp_1site_evolveNode()
%	Created by FS 23/02/2016

%% Perform TDVP on Node
% only use expvCustom due to complicated contractions!
fprintf('1');
para.sitej = 1; treeMPS.sitej = 1;
mps = tdvp_1site_evolveNode(mps,Vmat,para,op,results);			% nested function for now!
fprintf('-');



	function tdvp_1site_evolveNode()
	%% function tdvp_1site_evolveNode()
	%	evolves the Node = site 1 of a Tree-MPS
	%
	% Created 23/02/2016 by FS
	%
	t = para.tdvp.deltaT./2;
	
	if para.tdvp.imagT
		t = -1i*t;
	end
	
	d = size(treeMPS.mps{1});
	
	% skip time-evolution of Vmat for now, since only small dk!
	% add TDVP for Vmat if placing boson on node
	
	% update site 1 OBB terms h1j/h2jOBB
	if isempty(treeMPS.op.h1jOBB) || isempty(treeMPS.op.h1j)
		treeMPS.op = H_Eff(treeMPS, [], 'TR-A' , [], para);		% bring into OBB
	end
	
	%% Take matrix exponential
	% A (t+dt) = exp(-i ?? dt) * A(t)
	if para.tdvp.evolveSysTrotter == 0
		% old scheme, evolves system in one single step interacting with all chains at once
		[treeMPS.mps{1}, ~] = expvCustom(- 1i*t,'TREE-Hn1',...
								reshape(treeMPS.mps{1},[numel(treeMPS.mps{1}),1]), para, treeMPS.op);

		treeMPS.mps{1} = reshape(treeMPS.mps{1},d);
	else
		% newer scheme, trotterise in system-chain interactions to decrease complexity
		Atens = treeMPS.mps{1};
		dIn = d;
		for mc = 1:para.nChains
			%% create isometry into chain direction
			% rotate chain bonds
			para.currentChain = mc;
			if mc == 1
				[Atens,dOut] = tensShape(Atens, 'foldrotunfoldiso', 2, dIn);	% rot by 2 since leading singleton!
			else
				[Atens,dOut] = tensShape(Atens, 'foldrotunfoldiso', 1, dIn);	% rotates such that A: prod(D(1:NC without mc)) x (D(mc,1) * dk)
			end

			[Iso, A] = qr(Atens,0);			% Iso is isometry with all unused chains

			% TODO: comment if not testing!!
	% 		[Iso2,A2]= rrQR(Atens, floor(0.6*prod(dOut(end-1:end))),0);		% low rank QR approximation
	% 		fprintf('evolve: rank(A) = %d, dim(A,2) = %d, rrQR error = %g\n', rank(Atens), prod(dOut(end-1:end)),norm(Atens - Iso2*A2));
			
			treeMPS.tdvp.expvTol = para.tdvp.expvTol;
			treeMPS.tdvp.expvM   = para.tdvp.expvM;
			% evolve simplified mps matrix
			[A, err] = expvCustom(- 1i*t,'TREE-Hn1Trotter',...
				reshape(A,[numel(A),1]), para, op);
			A = reshape(A, prod(dOut(end-1:end)),[]);			% D*dk x D*dk
			% contract back together
			Atens = Iso * A;

			dIn = dOut;			% reset the new dimensions after rotation

		end

		mps{1} = reshape(Atens, d);
	end
	
	end

% function [A, Iso] = prepare_IsometryDkChain(Atens, para)
% 	%% Input: focused AS = mps{1} of star-MPS
% 	%  Splits Atens into Center A, containing bonds n_k, and D(para.currentChain)
% 	%	and the Isometry containing all other chain bonds
% 	%	see PEPS techniques for Trotter gates
% 	%  no truncation!
% 	% Output:
% 	%	A: D(NC,1)*dk x D(NC,1) x dk
% 	%	I:
% 	mc = para.currentChain;								% the chain to focus on
% 	d  = size(Atens);
%
% 	[AS]     = tensShape(Atens, 'unfold', mc+1, d);		% D(mc) x rest
%
% 	[AS, CA] = prepare_onesiteVmat(AS.', para);			% need for transpose since m x n, m > n input needed
%
% 	Atens	 = tensShape(AS.', 'fold', mc+1, d);
%
% end


end

function treeMPS = tdvp_1site_evolveChain(treeMPS,para)
%% function treeMPS = tdvp_1site_evolveChain(treeMPS,para)
%
%	evolves a single chain

%% Sweep r -> l:
for sitej = 1:paraChain.L
	fprintf('%g', sitej); paraChain.sitej = sitej;
	%% Update on-site Operators
	opChain = gen_sitej_op(opChain,paraChain,sitej,resultsChain.leftge);                     % take Site h1 & h2 Operators apply rescaling to Hleft, Hright, Opleft ...???

	%% Do the time-evolution of A and V with H(n)
	% this is symmetric for l->r and l<-r
	if ~paraChain.useVtens
		[mpsChain{sitej}, VmatChain{sitej}, paraChain, resultsChain, opChain, Hn] = tdvp_1site_evolveHn(mpsChain,VmatChain,paraChain,resultsChain,opChain,sitej);
	else
		paraChain.tdvp.expvCustomNow = 1;                                     % necessary setting!
		[mpsChain{sitej}, VmatChain{sitej}, paraChain, resultsChain, opChain, Hn] = tdvp_1site_evolveHnMC(mpsChain,VmatChain,paraChain,resultsChain,opChain,sitej);
	end
	% now: A and V are time-evolved, A is focused
	% OBBDim has increased by 50%. This must be truncated in the next
	% step again!
	% if sitej = L, then start lr sweep with decomposition of mps

	%% Take focus to next site, evolve with K(n)
	if sitej ~= paraChain.L
		fprintf('-');
		%% Left-normalize A and get Center matrix C(n,t+dt)_(rl,r)
		% expand/truncate BondDimensions here?
		[mpsChain{sitej}, Cn, paraChain,resultsChain] = prepare_onesite_truncate(mpsChain{sitej}, paraChain,sitej,resultsChain);

		%% Do the time-evolution of C
		% evolve non-site center between sitej and sitej+1
		% and focus A(n+1)
		[mpsChain, VmatChain, paraChain, resultsChain] = tdvp_1site_evolveKn(mpsChain,VmatChain,paraChain,resultsChain,opChain,sitej,Cn,Hn);
		clear('Hn','Cn');

		%% update Left Hamiltonian operators
		opChain = updateop(opChain,mpsChain,VmatChain,sitej,paraChain);

		if paraChain.useVmat
% 					truncateOBB(sitej);			% speedup by truncating within SVD from V to A ?
		end

	else % sitej = L
		%% Normalize with last SVD
		[mpsChain{sitej}, Cn] = prepare_onesite_truncate(mpsChain{sitej}, paraChain,sitej);
		% Cn = 1 approximately. can be thrown away -> mps normalized
		% finish with focus on A after time evolution

	end
end

% 		fprintf('\n');							%Debug!
% 		fprintf('%2g-',para.d_optnew);
% 		fprintf('\n');

%% Log vNE etc ?
% from prepare_onesite() and prepare_onesiteVmat():
% results.Amat_vNE  array
% results.Amat_sv   cell
% results.Vmat_vNE
% results.Vmat_sv
%% Output matrix dimensions if changed:
fprintf('\n');
if paraChain.tdvp.truncateExpandBonds
	fprintf('para.D:\t');
	out = strrep(mat2str(paraChain.D),';','\n');
	fprintf([out(2:end-1),'\n']);
	out = mat2str(cellfun(@(x) x(find(x,1,'last')),resultsChain.Amat_sv),2);
	fprintf([out(2:end-1),'\n']);
end
if paraChain.tdvp.expandOBB
	fprintf('para.d_opt:\t');
	out = strrep(mat2str(paraChain.d_opt),';','\n');
	fprintf([out(2:end-1),'\n']);
	out = mat2str(cellfun(@(x) x(end),resultsChain.Vmat_sv(end,~cellfun('isempty',resultsChain.Vmat_sv(end,:)))),2);
	fprintf([out(2:end-1),'\n']);
end

%% SWEEP l <- r:
paraChain.sweepto = 'l';
% now Hn = H(L)_(l'*r'*n',l*r*n)
for sitej = paraChain.L-1:-1:0
	fprintf('-'); paraChain.sitej = sitej;

	%% Right-normalize A(n+1) and get Center matrix C(n,t+dt)_(l,lr)
	% normalisation needed for updateop()!
	% Applies to bond n -> use sitej for result storage
	if paraChain.tdvp.truncateExpandBonds
		[mpsChain{sitej+1}, Cn, paraChain, resultsChain, sv, vNE] = prepare_onesite_truncate(mpsChain{sitej+1},paraChain,sitej+1, resultsChain);
	else
		[mpsChain{sitej+1}, Cn, paraChain, resultsChain, sv, vNE] = prepare_onesite(mpsChain{sitej+1},paraChain,sitej+1, resultsChain);
	end

	%% Do the time-evolution of C
	% evolve non-site center between sitej and sitej+1
	% and focus A(n)
	if sitej == 0
		[    Cn  , VmatChain, paraChain, resultsChain] = tdvp_1site_evolveKn(mpsChain,VmatChain,paraChain,resultsChain,opChain,sitej,Cn,Hn);
		% contract Cn with mps{1} after for loop
	else
		[mpsChain, VmatChain, paraChain, resultsChain] = tdvp_1site_evolveKn(mpsChain,VmatChain,paraChain,resultsChain,opChain,sitej,Cn,Hn);
	end
	clear('Hn');

	fprintf('%g', sitej);

	%% Expand / Truncate dk if needed!
	if para.tdvp.useDkExpand && sitej+1 >= 1
		% para is needed to call genh1h2term_onesite correctly
		[VmatChain{sitej+1}, para, paraChain, resultsChain, op, opChain] = truncateExpandDk(VmatChain{sitej+1},para, paraChain,resultsChain,op,opChain, sitej+1);
	end

	%% update right Hamiltonian operators
	% prepare operators in (sitej+1) for time-evolution on sitej
	paraChain.sitej = paraChain.sitej+1;					% needed for multi-chain reshape
	opChain         = updateop(opChain,mpsChain,VmatChain,sitej+1,paraChain);
	paraChain.sitej = paraChain.sitej-1;

	if paraChain.useVmat
% 				truncateOBB(sitej+1);
	end
	if sitej == 0
		continue;
	end
	%% Get on-site Operators and dimensions
	opChain = gen_sitej_op(opChain,paraChain,sitej,resultsChain.leftge);                     % take Site h1 & h2 Operators apply rescaling to Hleft, Hright, Opleft ...???

% 		shift = getObservable({'1siteshift'}, mps{sitej}, Vmat{sitej}, para); % calculate applicable shift
% 		fprintf('\n Shift: %s', mat2str(shift));
%
	%% Do the time-evolution of A and V
	% this is symmetric for l->r and l<-r
	if ~paraChain.useVtens
		[mpsChain{sitej}, VmatChain{sitej}, paraChain, resultsChain, opChain, Hn] = tdvp_1site_evolveHn(mpsChain,VmatChain,paraChain,resultsChain,opChain,sitej);
	else
		[mpsChain{sitej}, VmatChain{sitej}, paraChain, resultsChain, opChain, Hn] = tdvp_1site_evolveHnMC(mpsChain,VmatChain,paraChain,resultsChain,opChain,sitej);
	end
end

end


function L = getTreeLength(treeMPS)
%% function L = getTreeLength(treeMPS)
%	extracts the length of the treeMPS, based on the length of its children
%
if isfield(treeMPS,'L') && ~isempty(treeMPS.L)
	L = treeMPS.L;
	return
elseif treeMPS.height == 0
	if isfield(treeMPS,'mps')
		L = length(treeMPS.mps);
	else
		L = 0;					% is leaf without MPS
	end
	return
end

% first determine size(L)
h = treeMPS.height;			% ndims of L > 0
d = zeros(1,h);				% this will be size(L)
d(1) = treeMPS.degree+1;

if h == 1					% this is StarMPS with leaves only
	Lchild = arrayfun(@getTreeLength, treeMPS.child);
	L = [1;reshape(Lchild,d(1)-1,1)];
	return
end
% else height >= 2

if h == 2
	d(3) = 1;		% need singleton to allow comparisons later
end
Lchild = arrayfun(@getTreeLength, treeMPS.child,'UniformOutput',false);
% find out maximum dimensions per level to zero-pad Lchild
idx = find([treeMPS.child.height] == h-1);		% only pick the ones with maximum height
for ii = idx
	d(2:end) = max(d(2:end),size(Lchild{ii}));
end
% now d contains the maximum dimensions of L -> initialise
L = zeros(d);
L(1) = 1;					% length of nodes = 1 always for now!
for ii = 1:treeMPS.degree
	idx = arrayfun(@(x) 1:x,size(Lchild{ii}),'UniformOutput',false);		% create cell array containing index ranges of Lchild{ii}
	L(ii+1,idx{:}) = Lchild{ii};
end

end

