function tdvp_1site(mps,Vmat,para,results,op)
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
%   TODO:
%       - dk expansion / truncation

%% 0. Parameter settings:
if ~isfield(para,'tdvp')
    para.tdvp.tmax = 300;
        % For PPC:
        %   H defined in eV, h\bar left out
        %   -> real tmax = T * 6.58211928(15)×10^-16
    para.tdvp.deltaT = 10;                  % size of timeslice in units:
    para.tdvp.maxExpMDim = 10^2;            % For Lappy: 100, OE-PC: 80, pc52: 260; System dependent, use benchmark!
    para.tdvp.expvTol = 1e-9;               % error tolerance of expv(); default: 1e-7
    para.tdvp.expvM   = 50;                 % dim of Krylov subspace in expv(); default: 30
        % Sets threshold size for matrix exponential:
        %   if dim(A) < : use built-in expm(At)*v
        %   else        : use Expokit expv(t,A,v, expvTol, expvM)
        %   set maxExpMDim = 0 to only use expv()
    para.tdvp.rescaling = 0;                % turn on/off rescaling in TDVP
    para.rescaling = para.tdvp.rescaling;
end

% open file for writing results
% if storeMPS then para.tdvp.filename will contain tMPS
assert(isfield(para.tdvp,'filename'),'Need definition of para.tdvp.filename to save results!');
outFile = matfile(para.tdvp.filename,'Writable',true);

para.tdvp.smallFile = [para.tdvp.filename(1:end-4),'-small.mat'];
if isfield(para.tdvp,'storeMPS') && para.tdvp.storeMPS == 0
	para.tdvp.smallFile = para.tdvp.filename;
	clear('outFile');
end

resumeTDVP = 0;					% 0 = start new, 1 = resume, 2 = continue completed

try
	if exist('outFile','var')
		[startSize,~] = size(outFile,'tmps');					% = timeslice+1; = slice to start with now!
		if isprop(outFile,'currentSize')
			startSize = outFile.currentSize;
		end
		% resume incomplete TDVP
		resumeTDVP = 1;
		load(para.tdvp.smallFile, 'tresults');
	else
		%% find other way to determine resuming
		% perhaps by loading tresults and lastIdx?
		load(para.tdvp.smallFile, 'tresults');
		if exist('tresults','var') && isfield(tresults,'lastIdx')
			startSize = tresults.lastIdx;
			resumeTDVP = 1;
		else
			startSize = 1;				% fallback value. Perhaps change?
		end
	end
catch
	if size(mps,1) ~= 1 && size(Vmat,1) ~= 1
		% continue previously completed TDVP
		% since tmps and tVmat as input, but empty tdvp file.
		resumeTDVP = 2;
	end
end

% estimate the actual dimension in expm() from MPS dimensions.
% Not really needed now.
bondDim = [1 results.D{end} 1];
if para.useVmat
    para.tdvp.currentExpMDim = results.d_opt{end}.*bondDim(1:end-1).*bondDim(2:end);
else
    para.tdvp.currentExpMDim = results.dk{end}.*bondDim(1:end-1).*bondDim(2:end);
end

%% 1. sweep: condition the ground state MPS and prepare hamiltonian operator terms H(n) and K(n)
    % is already done in op.opstorage and op.hlrstorage
    %
    % also prepare other stuff
if resumeTDVP == 0
	if exist('outFile','var')
		% initialize tMPS and calc observables at the end
		fprintf('Initialize tMPS\n');
	% 	preallocation gives no benefit for cells?
	% 	outFile.tmps(length(para.tdvp.t),para.L) = {[]};		%
	% 	outFile.tVmat(length(para.tdvp.t),para.L) = {[]};
		outFile.tmps(1,para.L) = {[]};
		outFile.tVmat(1,para.L) = {[]};
		outFile.tmps(1,:) = mps;
		outFile.tVmat(1,:) = Vmat;
	else
		% calc tresults now to initialize it
		tresults = calTimeObservables(mps,Vmat,para);
	end
	if para.logging
		fprintf('Initialize results.tdvp\n');
		n = length(para.tdvp.t);
		results.tdvp.Vmat_sv	  = cell(n,para.L);
        results.tdvp.Vmat_sv(1,:) = results.Vmat_sv;
		results.tdvp.Amat_sv	  = cell(n,para.L-1);
        results.tdvp.Amat_sv(1,:) = results.Amat_sv;
		results.tdvp.d_opt		  = zeros(n,para.L);
        results.tdvp.d_opt(1,:)   = para.d_opt;
		results.tdvp.D   		  = zeros(n,para.L-1);
        results.tdvp.D(1,:)       = para.D;
		results.tdvp.dk 		  = zeros(n,para.L);
        results.tdvp.dk(1,:)      = para.dk;
		results.tdvp.expvTime     = [];
% 		results.tdvp.expError	  = zeros(length(para.tdvp.t)-1, (para.L-1)*4*2+3);
		results.tdvp.expError	  = zeros(length(para.tdvp.t)-1, 1);		% otherwise file gets too large
	end
	para.trustsite(end) = para.L;       % necessary to enable several truncation / expansion mechanisms
	% For fresh starts, define number of timeslices
	para.tdvp.slices = 1:(length(para.tdvp.t)-1);
	para.tdvp.calcTime = zeros(length(para.tdvp.t)-1,1);		% logs the elapsed time after each slice in h
elseif resumeTDVP == 1
	fprintf('Resuming aborted TDVP\n')
	para.tdvp.slices = startSize:(length(para.tdvp.t)-1);

	% need to regenerate hamiltonian operators?
% 	[para]		= SBM_genpara(para);					% chain parameters have not changed!
% 	[op,para]	= genh1h2term(para);					% Hamiltonian has not changed!
	[op]		= initstorage(mps, Vmat, op,para);
	% since op is backed up less often!
elseif resumeTDVP == 2
	%% If continue previously completed TDVP
	if exist('outFile','var');
		outFile.tmps = mps;
		outFile.tVmat = Vmat;
		outFile.tmps(length(para.tdvp.t),para.L) = {[]};
		outFile.tVmat(length(para.tdvp.t),para.L) = {[]};
	end
	startSize = size(mps,1);
	mps = mps(end,:);
	Vmat = Vmat(end,:);
	para.tdvp.slices = startSize:(length(para.tdvp.t)-1);
% not used anymore:
% 	if ~isprop(outFile,'mps')
% 	   outFile.mps = tmps(1,:);				% save the starting mps / Vmat, necessary for new files
% 	   outFile.Vmat = tVmat(1,:);
% 	end
	load(para.tdvp.fromFilename, 'tresults');
end

if (resumeTDVP ~= 0) && (size(results.tdvp.expError,2) > 1)
	% shrink down to 1 and make long enough
	expandBy = length(para.tdvp.t)-1 - size(results.tdvp.expError,1);
	results.tdvp.expError = [max(results.tdvp.expError,[],2); zeros(expandBy, 1)];	% take only the maximum error per sweep
end

if exist('tresults','var')
	save(para.tdvp.smallFile,'para','op','results','tresults','-v7.3');
else
	save(para.tdvp.smallFile,'para','op','results','-v7.3');
end

%% 2. time sweep

for timeslice = para.tdvp.slices
    para.sweepto = 'r';
	para.timeslice = timeslice;
	para.tdvp.deltaT = para.tdvp.t(timeslice+1)-para.tdvp.t(timeslice);
    fprintf('t = %g\n', para.tdvp.t(timeslice+1));
    % sweep l->r and time evolve each site
    for sitej = 1:para.L
        fprintf('%g', sitej);
        %% Update on-site Operators
        op = gen_sitej_op(op,para,sitej,results.leftge);                     % take Site h1 & h2 Operators apply rescaling to Hleft, Hright, Opleft ...???

        %% Do the time-evolution of A and V with H(n)
        % this is symmetric for l->r and l<-r
        [mps, Vmat, para, results, Hn] = tdvp_1site_evolveHn(mps,Vmat,para,results,op,sitej);

        % now: A and V are time-evolved, A is focused
        % OBBDim has increased by 50%. This must be truncated in the next
        % step again!
        % if sitej = L, then start lr sweep with decomposition of mps

        %% Take focus to next site, evolve with K(n)
        if sitej ~= para.L
            fprintf('-');
            %% Left-normalize A and get Center matrix C(n,t+dt)_(rl,r)
            % expand/truncate BondDimensions here?
            [mps{sitej}, Cn, para,results] = prepare_onesite_truncate(mps{sitej}, para,sitej,results);

            %% Do the time-evolution of C
            % evolve non-site center between sitej and sitej+1
            % and focus A(n+1)
            [mps, Vmat, para, results] = tdvp_1site_evolveKn(mps,Vmat,para,results,op,sitej,Cn,Hn);
            clear('Hn','Cn');

            if para.useVmat
                truncateOBB(sitej);
            end

            %% update Left Hamiltonian operators
            op = updateop(op,mps,Vmat,sitej,para);

        else % sitej = L
            %% Do Nothing.
            % finish with focus on A after time evolution
        end
    end

    %% Log vNE etc ?
    % from prepare_onesite() and prepare_onesiteVmat():
    % results.Amat_vNE  array
    % results.Amat_sv   cell
    % results.Vmat_vNE
    % results.Vmat_sv

    %% SWEEP l <- r:
    para.sweepto = 'l';
    % now Hn = H(L)_(l'*r'*n',l*r*n)
    for sitej = para.L-1:-1:1
        fprintf('-');

        %% Right-normalize A(n+1) and get Center matrix C(n,t+dt)_(l,lr)
        % normalisation needed for updateop()!
        % Applies to bond n -> use sitej for result storage
        [mps{sitej+1}, Cn, para] = prepare_onesite(mps{sitej+1},para,sitej+1);
        % if A SV needed: [mps{sitej+1}, Cn, para,results] = prepare_onesite(mps{sitej+1},para,sitej+1,results);

        %% Do the time-evolution of C
        % evolve non-site center between sitej and sitej+1
        % and focus A(n)
        [mps, Vmat, para, results] = tdvp_1site_evolveKn(mps,Vmat,para,results,op,sitej,Cn,Hn);
        clear('Hn','Cn');

        if para.useVmat
            truncateOBB(sitej+1);
        end

        fprintf('%g', sitej);

        %% update right Hamiltonian operators
        % prepare operators in (sitej+1) for time-evolution on sitej
        op = updateop(op,mps,Vmat,sitej+1,para);

        %% Get on-site Operators and dimensions
        op = gen_sitej_op(op,para,sitej,results.leftge);                     % take Site h1 & h2 Operators apply rescaling to Hleft, Hright, Opleft ...???

        %% Do the time-evolution of A and V
        % this is symmetric for l->r and l<-r
        [mps, Vmat, para, results, Hn] = tdvp_1site_evolveHn(mps,Vmat,para,results,op,sitej);

    end

    %% finished with both sweeps, focused on A(1) after time-evolution
	% report time estimates
    fprintf('\n');
	completePercent = round(timeslice./length(para.tdvp.t)*1000)./10;
	para.tdvp.calcTime(timeslice) = toc(para.tdvp.starttime)./3600;
	% TODO: might need optimisation for exp10 timesteps!
	n = 1:timeslice;
	n = n(para.tdvp.calcTime(n,1)>0);		% compensates for zeros due to resume of aborted TDVP
	hoursTotal = fnval(fnxtr(csaps([0,n],[0;para.tdvp.calcTime(n,1)])), length(para.tdvp.t)-1);
%  	hoursTotal = interp1([0,para.tdvp.calcTime(1:timeslice)],length(para.tdvp.t),'spline','extrap'); % bad extrapolation behaviour!
% 	hoursElapsed = toc(para.tdvp.starttime)./3600;
% 	hoursLeft = hoursElapsed./completePercent.*(100-completePercent);
	fprintf('Completed: %.3g%%, Time elapsed: %.2gh, Time left: %.2gh, Time total: %.2gh\n',...
		completePercent, para.tdvp.calcTime(timeslice),...
		hoursTotal-para.tdvp.calcTime(timeslice), hoursTotal);
    %% save tmps and tVmat and log parameters
	if exist('outFile','var')
		fprintf('Saving tMPS\n');
		outFile.tmps(timeslice+1,:)  = mps;					% writes to File
		outFile.tVmat(timeslice+1,:) = Vmat;				% writes to File
		outFile.currentSize = timeslice+1;
	end
	if para.logging
        results.tdvp.Vmat_sv(timeslice+1,:) = results.Vmat_sv;
        results.tdvp.Amat_sv(timeslice+1,:) = results.Amat_sv;
        results.tdvp.d_opt(timeslice+1,:)   = para.d_opt;
        results.tdvp.D(timeslice+1,:)       = para.D;
        results.tdvp.dk(timeslice+1,:)      = para.dk;
	end
	if ~exist('tresults','var')
		% calculate everything up to now!
		if exist('outFile','var')
			tresults = calTimeObservables(outFile.tmps(1:outFile.currentSize,:),outFile.tVmat(1:outFile.currentSize,:),para);
		else
			tresults = calTimeObservables(mps,Vmat,para);
		end
	elseif isfield(tresults,'lastIdx') && tresults.lastIdx == timeslice
		% only calculate the current slice
		tresults = calTimeObservables(mps,Vmat,para,tresults);
	end
	%% Output matrix dimensions:
	fprintf('para.D:\n');
	fprintf('%2g-',para.D);
	fprintf('\npara.d_opt:\n');
	fprintf('%2g-',para.d_opt);
	fprintf('\n');
	%% Additional saving
	if mod(timeslice,para.tdvp.saveInterval) == 0 || timeslice ==  para.tdvp.slices(end)
		%% backup results, op and para less often than tmps and tVmat
		results.tdvp.time = toc(para.tdvp.starttime);
		fprintf('Saving results.. ')
		copyfile(para.tdvp.smallFile,[para.tdvp.smallFile(1:end-4),'.bak']);
		save(para.tdvp.smallFile,'para','op','results','tresults');
		if isfield(para.tdvp,'storeMPS') && para.tdvp.storeMPS == 0
			fprintf('and MPS..')
			save(para.tdvp.smallFile,'mps','Vmat','-append');
			save([para.tdvp.filename(1:end-4),'-small.mat'],'para','tresults');		% still need smaller file since results grows very much!
		end
		fprintf('\n');
	end
end

    function truncateOBB(sitej)
        %% Truncate OBB of mps{sitej}
        % code taken from adjustdopt.m
        % only apply if also expanded OBB before
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

end