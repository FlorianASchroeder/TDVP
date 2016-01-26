function tresults = calTimeObservables(tmps,tVmat,para,varargin)
%% function tresults = calTimeObservables(tmps,tVmat,para,varargin)
%	creates tresults and computes Observables for each timeslice in tmps
%
%	tresults = calTimeObservables(tmps,tVmat,para,tresults)
%		Appends the slices in tmps to the end of tresults.
%			e.g.: tresults ranges in t = [0:3]; tmps has slices t=[4:6]
%				  then returned tresults will have t = [0:6]
%
% Modified:
%	FS 23/07/15: - changed nx to n for Multi-Chain compatibility!
%	FS 24/07/15: - added switch to allow observable selection (especially current)
%   FS 18/08/15: - changed to single precision for smaller tresults file!

	if isfield(para.tdvp,'extractObsInterval')
		% only works with equidistant steps and single tmps slices
		if mod(para.tdvp.tmax, para.tdvp.extractObsInterval) == 0 && (para.tdvp.extractObsInterval >= para.tdvp.deltaT)
			totalN = round(para.tdvp.tmax/para.tdvp.extractObsInterval) +1;
		else
			error('VMPS:calTimeObservables:WrongExtractObsInterval','Need to define extractObsInterval so that mod(tmax,interval)=0!');
		end
	else
		totalN = size(para.tdvp.t,2);
	end
	if nargin > 3
		% continue tresults
		tresults = varargin{1};
		assert(isfield(tresults,'t'),'4th argument has to be tresults');
		if ~isfield(tresults,'lastIdx')			% for compatibility with old files
			tresults.lastIdx = size(tresults.n,1);
		end
		missingN = totalN - size(tresults.n,1);
		if missingN > 0
			% pre-allocate memory
			tresults.n(totalN,end,max(para.nChains,para.nEnvironments)) = single(0);
			tresults.t(1,totalN)	= single(0);
		end
	else
		tresults.lastIdx = 0; missingN = 0;
		fprintf('Calculate Observables:\n');
		tresults.n  = single(zeros(totalN,para.L,max(para.nChains,para.nEnvironments)));
		tresults.t  = single(zeros(1,totalN));
		para.timeslice = 0;						% needed for extractObsInterval first run
	end

	slices = (1:size(tmps,1))+tresults.lastIdx;		% in tresults file

	if isfield(para.tdvp,'extractObsInterval')
		if mod(para.tdvp.t(1,para.timeslice+1),para.tdvp.extractObsInterval) ~= 0
			return;		% skip this timeslice!
		end
	end

	for i = slices
		if length(slices) > 1
			fprintf('%g-',i);
		else
			fprintf('<O> for slice %d\n',i);
		end
		j = i-tresults.lastIdx;			% shifted index for tmps, tVmat

        %% General Observables
        % 1. Chain Occupation
		tresults.n(i,:,:) = single(getObservable({'occupation'},tmps(j,:),tVmat(j,:),para));     % (L x nChain)
		tresults.t(i)     = single(para.tdvp.t(1,para.timeslice+1));

		% 2. Boson chain observables
		if strContains(para.tdvp.Observables,'.j.','.sn.')
			% save here, to reuse later!
			AnAm = getObservable({'bath2correlators'}, tmps(j,:),tVmat(j,:),para);
		end

		% Chain polaron
		if strContains(para.tdvp.Observables,'.x.','.x2.')		% diabatic or adiabatic states
			% 2*Re<a^+>
			if strfind(para.tdvp.Observables,'.x.')
				chainX = 2*real(getObservable({'bath1correlators'}, tmps(j,:),tVmat(j,:),para));				% L x nStates x nChains
			elseif strfind(para.tdvp.Observables,'.x2.')
				chainX = 2*real(getObservable({'bath1correlators','adiabatic'}, tmps(j,:),tVmat(j,:),para));	% L x nStates x nChains
			end
			[~, nStates, nChains] = size(chainX);
			if ~isfield(tresults,'x') || missingN > 0
				tresults.x(totalN,para.L,nStates, nChains) = single(0);
			end
			tresults.x(i,:,:,:) = chainX;
		end

		if strfind(para.tdvp.Observables,'.j.')
			%% Calculate current along entire chain
			if ~isfield(tresults,'j') || missingN > 0
				tresults.j(totalN,para.L,nChains) = single(0);
			end
			if ~exist('AnAm','var')
				tresults.j(i,:,:) = single(getObservable({'current'},tmps(j,:),tVmat(j,:),para));
			else
				tresults.j(i,:,:) = single(getObservable({'current',AnAm},tmps(j,:),tVmat(j,:),para));
			end
		end

		% 3. Star Observables
		if isfield(para.tdvp,'extractStarInterval') && strContains(para.tdvp.Observables,'.sn.','.sx.','.sx2.')
			Nslice = round(para.tdvp.extractStarInterval / para.tdvp.extractObsInterval);		% how often to extract Star Observables
			if mod(i-1,Nslice) == 0
				pos = ceil(i/Nslice);

				if strfind(para.tdvp.Observables,'.sn.')
					if ~exist('AnAm','var')
						occ	= getObservable({'staroccupation'}     ,tmps(j,:),tVmat(j,:),para);		% (1+1) x k x nc
					else
						occ = getObservable({'staroccupation',AnAm},tmps(j,:),tVmat(j,:),para);		% (1+1) x k x nc
					end
					starOmega = squeeze(single(occ(1,:,:)));				% get rid of leading singleton
				end
				
				if strContains(para.tdvp.Observables,'.sx.','.sx2.')
					if strfind(para.tdvp.Observables,'.sx.')
						polaron = getObservable({'starpolaron'},tmps(j,:),tVmat(j,:),para);				% (1+2) x k x nc, diabatic states
					elseif strfind(para.tdvp.Observables,'.sx2.')
						% not state projecting, but selecting single dominating states across the first bond!
						% kind of similar to diabatic states picture!
						polaron = getObservable({'starpolaron','adiabatic'},tmps(j,:),tVmat(j,:),para);	% (1+2) x k x nc

					end
					starOmega = squeeze(single(polaron(1,:,:)));			% get rid of leading singleton
				end
				
				if ~isfield(tresults, 'star')
					% initialise storage if first sweep
					nElements = para.tdvp.tmax/para.tdvp.extractStarInterval +1;
					if exist('occ','var'),     tresults.star.n = single(zeros(nElements,size(occ,2),para.nChains)); end
					if exist('polaron','var'), tresults.star.x = single(zeros(nElements,size(polaron,2),2,para.nChains)); end
					tresults.star.omega = starOmega;
					tresults.star.t     = single(zeros(1,nElements));
				end

				if strfind(para.tdvp.Observables,'.sn.')
					tresults.star.n(pos,:,:) = single(occ(2,:,:));
				end
				if strContains(para.tdvp.Observables,'.sx.','.sx2.')
					for kk = 1:size(polaron,1)-1
						tresults.star.x(pos,:,kk,:) = single(polaron(kk+1,:,:));	% SBM: 1: up-proj, 2: down-proj
					end
				end
				tresults.star.t(pos)   = single(para.tdvp.t(1,para.timeslice+1));
			end
		end

		if ~isempty(strfind(para.model, 'SpinBoson'))
            %% Observables for SBM
            % 1. Spin Observables
			if ~isfield(tresults,'spin')
                tresults.spin.sx = single(zeros(totalN,1));
                tresults.spin.sy = single(zeros(totalN,1));
                tresults.spin.sz = single(zeros(totalN,1));
                tresults.spin.visibility = single(zeros(totalN,1));
			elseif missingN > 0
				tresults.spin.sx = single([tresults.spin.sx; zeros(missingN,1)]);
				tresults.spin.sy = single([tresults.spin.sy; zeros(missingN,1)]);
				tresults.spin.sz = single([tresults.spin.sz; zeros(missingN,1)]);
				tresults.spin.visibility = single([tresults.spin.visibility; zeros(missingN,1)]);
			end
            temp = getObservable({'spin'},tmps(j,:),tVmat(j,:),para);
            tresults.spin.sx(i) = single(temp.sx);
            tresults.spin.sy(i) = single(temp.sy);
            tresults.spin.sz(i) = single(temp.sz);
            tresults.spin.visibility(i) = single(sqrt(temp.sx^2+temp.sy^2));
		end

		% 1. Density matrix
		if strContains(para.tdvp.Observables,'.dm.','.dm2.')
			if ~isfield(tresults,'rho')
				tresults.rho = single(zeros(totalN,para.dk(1,1),para.dk(1,1)));
			elseif missingN > 0
				tresults.rho = single([tresults.rho; zeros(missingN,para.dk(1,1),para.dk(1,1))]);
			end
			tresults.rho(i,:,:,1) = single(getObservable({'rdm',1},tmps(j,:),tVmat(j,:),para));
			if strContains(para.tdvp.Observables,'.dm2.')
				% only for 2-lvl system for now; only calculates largest bond state.
				tresults.rho(i,:,:,2) = single(getObservable({'rdm_adiabatic',1,1},tmps(j,:),tVmat(j,:),para));  %{'rdm_adiabatic',sitej,state}
				tresults.rho(i,:,:,3) = single(getObservable({'rdm_adiabatic',1,2},tmps(j,:),tVmat(j,:),para));  %{'rdm_adiabatic',sitej,state}
			end
		end

		if strcmp(para.model, 'MLSBM') || ~isempty(strfind(para.model,'DPMES'))
            %% Observables for MLSBM
            % 2. PPC Wavefunction
			%    only if not extracting Density Matrix
			if ~strfind(para.tdvp.Observables,'.dm.')
				if ~isfield(tresults,'PPCWavefunction')
					tresults.PPCWavefunction = single(zeros(totalN,para.dk(1,1)));
				elseif missingN > 0
					tresults.PPCWavefunction = single([tresults.PPCWavefunction; zeros(missingN,para.dk(1,1))]);
				end
				tresults.PPCWavefunction(i,:) = single(diag(getObservable({'rdm',1},tmps(j,:),tVmat(j,:),para)));
			end

            % 3. Participation on ring
			if ~isfield(tresults,'participation')
                tresults.participation = single(zeros(totalN,1));
			elseif missingN > 0
				tresults.participation = single([tresults.participation; zeros(missingN,1)]);
			end
            tresults.participation(i) = single(getObservable({'participation'},tmps(j,:),tVmat(j,:),para));

			% 4. Hs + Hi
			if ~isfield(tresults,'hshi')
				tmp = single(getObservable({'hshi'},tmps(j,:),tVmat(j,:),para));
                tresults.hshi = single(zeros(totalN,length(tmp)));
				tresults.hshi(i,:) = tmp;
			elseif missingN > 0
				tresults.hshi = single([tresults.hshi; zeros(missingN,1)]);
			else
				tresults.hshi(i,:) = single(getObservable({'hshi'},tmps(j,:),tVmat(j,:),para));
			end

		end
		
		if strContains(para.tdvp.Observables,'.sp.')					% sp for state projection
			if ~isfield(tresults,'stateProjection')
				tresults.stateProjection = single(zeros(totalN,1));
			elseif missingN > 0
				tresults.stateProjection(totalN,1) = 0;			% does preallocation
			end
			tresults.stateProjection(i,1) = single(getObservable({'stateproject',para.InitialState,1},tmps(j,:),tVmat(j,:),para));	% project onto |IS>|0>, IS = initial state
		end
		
		if strContains(para.tdvp.Observables,'.ss.')					% ss for system state
			if ~isfield(tresults,'system') || ~isfield(tresults.system,'state')
				tresults.system.state = single(zeros(totalN,para.dk(1),para.dk(1)));		% t x dk x D (adiabatic)
			elseif missingN > 0
				tresults.system.state(totalN,para.dk(1),para.dk(1)) = 0;			% does preallocation
			end
			tresults.system.state(i,:,:) = single(getObservable({'state',1},tmps(j,:),tVmat(j,:),para));
		end
		
		if strContains(para.tdvp.Observables,'.ses.')					% ses for system-environment state
			if ~isfield(tresults,'mps')
				tresults.mps = cell(totalN,2);		% t x sites
				tresults.Vmat = cell(totalN,2);		% t x sites
			elseif missingN > 0
				tresults.mps(totalN,2) = {};			% does preallocation
				tresults.Vmat(totalN,2) = {};
			end
			out = getObservable({'sys-env-state'},tmps(j,:),tVmat(j,:),para);		% get mps([1,2]) and Vmat([1,2])
			tresults.mps(i,:) = out.mps;
			tresults.Vmat(i,:) = out.Vmat;
		end
		
		if strContains(para.model, 'SpinBosonTTM', 'UniformBosonTTM')
			%% extract transfer tensor
			% only use for single-slice tMPS due to iterative procedure
			if length(slices) > 1, return; end;
			rdm = getObservable({'rdm',[1 2]},tmps(j,:),tVmat(j,:),para);
			% Ortho Normal Operator Basis in dxd
			d = para.dk(1,2);
			ONOB = eye(d^2); ONOB = reshape(ONOB,[d,d,d^2]);
			EAm = zeros(d,d,d^2);
			if d <= 4
				for k = 1:d^2
					% slow ncon
					EAm(:,:,k) = ncon({rdm, squeeze(ONOB(:,:,k))'},...
									  {[-1,2,-2,1], [1,2]})*d;			% apply Op, contract / trace; perhaps *d
					% EAm_ijk  = d* rdm_imjn * ONOB*_mnk
				end
			else
				EAm = d*contracttensors(rdm, 4, [2 4], conj(ONOB), 3, [1 2]);	% EAm_ijk  = d* rdm_imjn * ONOB*_mnk
			end
			Epsilon = reshape(EAm,[d^2,d^2]);	% E_(n~',n~),(n',n)
			T = Epsilon;
% 			for k = 3:i
% 				T = T - tresults.TTM.T(:,:,i+1-k)*tresults.TTM.Epsilon(:,:,k-1);
% 			end
			for k = 2:i-1
				if ~para.tdvp.splitTTM
					T = T - tresults.TTM.T(:,:,i-k)*tresults.TTM.Epsilon(:,:,k);			% just did k -> k+1
				else
% 					Test = T - tresults.TTM.T(:,:,i-k)*tresults.TTM.Epsilon(:,:,k);			% just did k -> k+1
					% looks complicated, but should have better scaling O(d^5*r) than above O(d^6)! (for rank r approx.)
% 					tempE = reshape(tresults.TTM.Epsilon(:,:,k),d,d,d^2);					% E_(n',n,j)
% 					tempE = contracttensors(conj(tresults.TTM.TV{i-k}),3,2, tempE,3,2);		% E_(n~,r,n',j) = V*_(n~,n,r) * E_(n',n,j)
% 					tempE = contracttensors(tresults.TTM.TD{i-k},2,2,tempE,4,2);			% E_(r',n~,n',j) = D_(r',r) * E_(n~,r,n',j)
% 					tempE = contracttensors(tresults.TTM.TV{i-k}, 3, [2,3], tempE,4,[3,1]);	% E_(n~',n~,j) = V_(n~',n',r') * E_(r',n~,n',j)
% 					T = T - reshape(tempE,[d^2,d^2]);
% 					
					% MUCH faster!! perhaps uses symmetries! As fast as without splitting
					tempT = tresults.TTM.TV{i-k}*sparse(diag(tresults.TTM.TD{i-k}))*tresults.TTM.TV{i-k}';% T_(n~',n'),(n~,n)
					tempT = reshape(tempT,d,d,d,d);
					tempT = permute(tempT,[1,3,2,4]);
					tempT = reshape(tempT,d^2,d^2);
					T = T - tempT * tresults.TTM.Epsilon(:,:,k);
				end
			end
			% following vectorisation is slower than for-loop! due to inverse ordering??
% 			if i >= 3
% 				T = T - contracttensors(tresults.TTM.T(:,:,i-(2:i-1)),3,[2 3], tresults.TTM.Epsilon(:,:,2:(i-1)),3,[1 3]);
% 			end
			tresults.TTM.Epsilon(:,:,i) = Epsilon;
			if i > 1
				if ~para.tdvp.splitTTM
					tresults.TTM.Tnorm(i-1) = single(norm(T));
					tresults.TTM.T(:,:,i-1) = T;
				else
% 					tresults.TTM.T(:,:,i-1) = T;
					% reshape T to allow low-rank approximation
					T = reshape(T,[d,d,d,d]);		% T_(n~',n~,n',n)
					T = permute(T,[1,3,2,4]);		% T_(n~',n',n~,n)
					T = reshape(T,[d^2,d^2]);		% T_(n~',n'),(n~,n)
					
					[V,D] = eig(T);					% do eig to obtain proper self-adjointness
					D = real(diag(D));				% T self-adjoint -> D is real!
% 					plot(abs(D),'Displayname',num2str(i)); hold all; set(gca,'Yscale','log'); drawnow
					keepdims = abs(D)>1e-14;		% safe threshold?
					V = V(:,keepdims);
					D = diag(D(keepdims));			% now: norm(T-V*D*V') < 1e-14
					% obtain T now as V*D*V'; V_(n~',n'),i
% 					tresults.TTM.TV{i-1} = reshape(V,d,d,sum(keepdims));
					tresults.TTM.TV{i-1} = V;		% d^2 x r
					tresults.TTM.TD{i-1} = diag(D);	% store as vector
					tresults.TTM.Tnorm(i-1) = norm(D);
				end
				fprintf('\n|TTM|/dt^2: %g\n',tresults.TTM.Tnorm(i-1)/para.tdvp.deltaT^2);
			end
		end
		missingN = 0;					% shall only be used in first loop!
	end
	tresults.lastIdx = i;
    fprintf('\n');
end

function out = strContains(str, varargin)
out = false;
for ii = 1:length(varargin)
	out = out || ~isempty(strfind(str,varargin{ii}));
end

end