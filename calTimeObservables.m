function tresults = calTimeObservables(tmps,tVmat,para,varargin)
%% function tresults = calTimeObservables(tmps,tVmat,para,varargin)
%	creates tresults and computes Observables for each timeslice in tmps
%
%	tresults = calTimeObservables(tmps,tVmat,para,tresults)
%		Appends the slices in tmps to the end of tresults.
%			e.g.: tresults ranges in t = [0:3]; tmps has slices t=[4:6]
%				  then returned tresults will have t = [0:6]

	if isfield(para.tdvp,'extractObsInterval')
		% only works with equidistant steps and single tmps slices
		if mod(para.tdvp.tmax, para.tdvp.extractObsInterval) == 0 && (para.tdvp.extractObsInterval >= para.tdvp.deltaT)
			totalN = para.tdvp.tmax/para.tdvp.extractObsInterval +1;
		else
			error('Need to define extractObsInterval so that mod(tmax,interval)=0!');
		end
	else
		totalN = size(para.tdvp.t,2);
	end
	if nargin > 3
		% continue tresults
		tresults = varargin{1};
		assert(isfield(tresults,'nx'),'4th argument has to be tresults');
		if ~isfield(tresults,'lastIdx')			% for compatibility with old files
			tresults.lastIdx = size(tresults.nx,1);
		end
		missingN = totalN - size(tresults.nx,1);
		if missingN > 0
			tresults.nx(totalN,end) = 0;
			tresults.t(1,totalN)	= 0;
		end
	else
		tresults.lastIdx = 0; missingN = 0;
		fprintf('Calculate Observables:\n');
	    tresults.nx = zeros(totalN,para.L);
		tresults.t  = zeros(1,totalN);
		para.timeslice = 0;						% needed for extractObsInterval first run
	end

	slices = (1:size(tmps,1))+tresults.lastIdx;		% in tresults file

	if isfield(para.tdvp,'extractObsInterval')
		if mod(para.tdvp.t(1,para.timeslice+1),para.tdvp.extractObsInterval) ~= 0
			return;		% skip this timeslice!
		end
	end

	for i = slices
        fprintf('%g-',i);
		j = i-tresults.lastIdx;			% shifted index for tmps, tVmat

        %% General Observables
        % 1. Chain Occupation
        tresults.nx(i,:) = getObservable({'occupation'},tmps(j,:),tVmat(j,:),para);
		tresults.t(i)    = para.tdvp.t(1,para.timeslice+1);

		% 2. Star Occupation
		if isfield(para.tdvp,'extractStarInterval')
			if mod(para.tdvp.t(1,para.timeslice+1),para.tdvp.extractStarInterval) == 0
				pos = round(para.tdvp.t(1,para.timeslice+1)/para.tdvp.extractStarInterval) +1;

				occ		= getObservable({'staroccupation'},tmps(j,:),tVmat(j,:),para);		% 2 x k
				polaron = getObservable({'starpolaron'},tmps(j,:),tVmat(j,:),para);			% 2 x k

				if ~isfield(tresults, 'star')
					% initialise storage if first sweep
					nElements = para.tdvp.tmax/para.tdvp.extractStarInterval +1;
					tresults.star.n	    = zeros(nElements,length(occ));
					tresults.star.x	    = zeros(nElements,length(polaron));
					tresults.star.omega = occ(1,:);
					tresults.star.t     = zeros(1,nElements);
				end

				tresults.star.n(pos,:) = occ(2,:);
				tresults.star.x(pos,:) = polaron(2,:);
				tresults.star.t(pos)   = para.tdvp.t(1,para.timeslice+1);
			end
		end

		if strcmp(para.model, 'SpinBoson')
            %% Observables for SBM
            % 1. Spin Observables
			if ~isfield(tresults,'spin')
                tresults.spin.sx = zeros(totalN,1);
                tresults.spin.sy = zeros(totalN,1);
                tresults.spin.sz = zeros(totalN,1);
                tresults.spin.visibility = zeros(totalN,1);
			elseif missingN > 0
				tresults.spin.sx = [tresults.spin.sx; zeros(missingN,1)];
				tresults.spin.sy = [tresults.spin.sy; zeros(missingN,1)];
				tresults.spin.sz = [tresults.spin.sz; zeros(missingN,1)];
				tresults.spin.visibility = [tresults.spin.visibility; zeros(missingN,1)];
			end
            temp = getObservable({'spin'},tmps(j,:),tVmat(j,:),para);
            tresults.spin.sx(i) = temp.sx;
            tresults.spin.sy(i) = temp.sy;
            tresults.spin.sz(i) = temp.sz;
            tresults.spin.visibility(i) = sqrt(temp.sx^2+temp.sy^2);
		end

		if strcmp(para.model, 'MLSBM')
            %% Observables for MLSBM
            % 2. PPC Wavefunction
			if ~isfield(tresults,'PPCWavefunction')
                tresults.PPCWavefunction = zeros(totalN,16);
			elseif missingN > 0
				tresults.PPCWavefunction = [tresults.PPCWavefunction; zeros(missingN,16)];
			end
            tresults.PPCWavefunction(i,:) = diag(getObservable({'rdm',1},tmps(j,:),tVmat(j,:),para));

            % 3. Participation on ring
			if ~isfield(tresults,'participation')
                tresults.participation = zeros(totalN,1);
			elseif missingN > 0
				tresults.participation = [tresults.participation; zeros(missingN,1)];
			end
            tresults.participation(i) = getObservable({'participation'},tmps(j,:),tVmat(j,:),para);
		end
		missingN = 0;					% shall only be used in first loop!
	end
	tresults.lastIdx = i;
    fprintf('\n');
end