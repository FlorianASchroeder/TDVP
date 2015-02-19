function tresults = calTimeObservables(tmps,tVmat,para,varargin)
%% function tresults = calTimeObservables(tmps,tVmat,para,varargin)
%	creates tresults and computes Observables for each timeslice in tmps
%
%	tresults = calTimeObservables(tmps,tVmat,para,tresults)
%		Appends the slices in tmps to the end of tresults.
%			e.g.: tresults ranges in t = [0:3]; tmps has slices t=[4:6]
%				  then returned tresults will have t = [0:6]

	totalN = size(para.tdvp.t,2);
	if nargin > 3
		tresults = varargin{1};
		assert(isfield(tresults,'nx'),'4th argument has to be tresults');
		if ~isfield(tresults,'lastIdx')			% for compatibility with old files
			tresults.lastIdx = size(tresults.nx,1);
		end
		missingN = totalN - size(tresults.nx,1);
		tresults.nx = [tresults.nx; zeros(missingN,para.L)];
	else
		tresults.lastIdx = 0; missingN = 0;
		fprintf('Calculate Observables:\n');
	    tresults.nx = zeros(totalN,para.L);
	end

	slices = (1:size(tmps,1))+tresults.lastIdx;

    for i = slices
        fprintf('%g-',i);
		j = i-tresults.lastIdx;			% shifted index for tmps, tVmat

        %% General Observables
        % 1. Chain Occupation
        tresults.nx(i,:) = getObservable({'occupation'},tmps(j,:),tVmat(j,:),para);

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