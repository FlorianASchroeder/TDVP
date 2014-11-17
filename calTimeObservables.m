function tresults = calTimeObservables(tmps,tVmat,para)
    fprintf('Calculate Observables:\n');

    tresults.nx = zeros(size(tmps,1),para.L);

    for i = 1:size(tmps,1)
        fprintf('%g-',i);

        %% General Observables
        % 1. Chain Occupation
        tresults.nx(i,:) = getObservable({'occupation'},tmps(i,:),tVmat(i,:),para);

        if strcmp(para.model, 'SpinBoson')
            %% Observables for SBM
            % 1. Spin Observables
            if ~isfield(tresults,'spin')
                tresults.spin.sx = zeros(size(tmps,1),1);
                tresults.spin.sy = zeros(size(tmps,1),1);
                tresults.spin.sz = zeros(size(tmps,1),1);
                tresults.spin.visibility = zeros(size(tmps,1),1);
            end
            temp = getObservable({'spin'},tmps(i,:),tVmat(i,:),para);
            tresults.spin.sx(i) = temp.sx;
            tresults.spin.sy(i) = temp.sy;
            tresults.spin.sz(i) = temp.sz;
            tresults.spin.visibility(i) = sqrt(temp.sx^2+temp.sy^2);
        end


        if strcmp(para.model, 'MLSBM')
            %% Observables for MLSBM
            % 2. PPC Wavefunction
            if ~isfield(tresults,'PPCWavefunction')
                tresults.PPCWavefunction = zeros(size(tmps,1),16);
            end
            tresults.PPCWavefunction(i,:) = diag(getObservable({'rdm',1},tmps(i,:),tVmat(i,:),para));

            % 3. Participation on ring
            if ~isfield(tresults,'participation')
                tresults.participation = zeros(size(tmps,1),1);
            end
            tresults.participation(i) = getObservable({'participation'},tmps(i,:),tVmat(i,:),para);
        end
    end
    fprintf('\n');
end