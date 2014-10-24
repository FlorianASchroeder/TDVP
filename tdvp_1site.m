function [mps, Vmat, para, tmps, tVmat] = tdvp_1site(mps,Vmat,para,results,op)
%% VMPS Time-dependent variational principle
%   implementation following Haegeman et al. 2014, arXiv:1408.5056
%   using Expokit for expv()
%
%   (only use without OBB (Vmat) at the moment!)
%
%   Created by Florian Schroeder 07/10/2014
%
%

%% 0. Parameter settings:

para.tdvp.tmax = 300;
    % For PPC:
    %   H defined in eV, h\bar left out
    %   -> real tmax = T * 6.58211928(15)×10^-16
para.tdvp.deltaT = 10;                 % size of timeslice in units:
para.tdvp.maxExpMDim = 10^2;            % maximum allowed a_n*a_(n+1)*d_k.
para.tdvp.expvTol = 1e-9;               % error tolerance of expv(); default: 1e-7
para.tdvp.expvM   = 50;                 % dim of Krylov subspace in expv(); default: 30
    % Sets threshold size for matrix exponential:
    %   if dim(A) < : use built-in expm(At)*v
    %   else        : use Expokit expv(t,A,v, expvTol, expvM)
    %   set maxExpMDim = 0 to only use expv()
para.tdvp.rescaling = 0;                % turn on/off rescaling in TDVP
para.rescaling = para.tdvp.rescaling;

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

    % initial time values
    tmps(1, :) = mps;
    tVmat(1,:) = Vmat;
%% 2. time sweep

for timeslice = 1:(para.tdvp.tmax/para.tdvp.deltaT)
    para.sweepto = 'r';
    fprintf('t = %g\n', timeslice * para.tdvp.deltaT);
    % sweep l->r and time evolve each site
    for sitej = 1:para.L
        fprintf('%g', sitej);
        %% Update on-site Operators
        op = gen_sitej_op(op,para,sitej,results.leftge);                     % take Site h1 & h2 Operators apply rescaling to Hleft, Hright, Opleft ...???

        %% Do the time-evolution of A and V with H(n)
        % this is symmetric for l->r and l<-r
        [mps, Vmat, para, results, Hn] = tdvp_1site_evolveHn(mps,Vmat,para,results,op,sitej);

        % now: A and V are time-evolved, A is focused
        % if sitej = L, then start lr sweep with decomposition of mps

        %% Take focus to next site, evolve with K(n)
        if sitej ~= para.L
            fprintf('-');
            %% Left-normalize A and get Center matrix C(n,t+dt)_(rl,r)
            [mps{sitej}, Cn, para,results] = prepare_onesite(mps{sitej}, para,sitej,results);

            %% Do the time-evolution of C
            % evolve non-site center between sitej and sitej+1
            % and focus A(n)
            [mps, Vmat] = tdvp_1site_evolveKn(mps,Vmat,para,results,op,sitej,Cn,Hn);
            clear('Hn','Cn');


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
        [mps{sitej+1}, Cn, para,results] = prepare_onesite(mps{sitej+1},para,sitej,results);

        %% Do the time-evolution of C
        % evolve non-site center between sitej and sitej+1
        % and focus A(n)
        [mps, Vmat] = tdvp_1site_evolveKn(mps,Vmat,para,results,op,sitej,Cn,Hn);
        clear('Hn','Cn');

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

    % finished with both sweeps, focused on A(1) after time-evolution
    fprintf('\n');
    %% calculate time-dependent expectation values
    tmps(timeslice+1, :) = mps;
    tVmat(timeslice+1,:) = Vmat;
end


end