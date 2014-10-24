%% Do Time-Evolution with 1-site TDVP
[mps1, Vmat1, para, tmps, tVmat] = tdvp_1site(mps,Vmat,para,results,op);

%% calculate observables:

tresults.nx = zeros(size(tmps,1),para.L);
tresults.participation = zeros(size(tmps,1),1);
tresults.PPCWavefunction = zeros(size(tmps,1),16);
for i = 1:size(tmps,1)
    fprintf('%g-',i);
    % 1. Occupation
    tresults.nx(i,:) = getObservable({'occupation'},tmps(i,:),tVmat(i,:),para);

    % 2. Participation on ring
    tresults.participation(i) = getObservable({'participation'},tmps(i,:),tVmat(i,:),para);

    % 3. PPC Wavefunction
    if strcmp(para.model, 'MLSBM')
        tresults.PPCWavefunction(i,:) = diag(getObservable({'rdm',1},tmps(i,:),tVmat(i,:),para));
    end
end

%%
