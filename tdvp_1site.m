function [mps, Vmat] = tdvp_1site(mps,Vmat,para,results,op)
%% VMPS Time-dependent variational principle
%   implementation following Haegeman et al. 2014, arXiv:1408.5056
%
%   (only use without OBB (Vmat) at the moment!)
%
%   Created by Florian Schroeder 07/10/2014
%
%

%% 0. Parameter settings:

para.tdvp.tmax = 1;
para.tdvp.deltaT = 0.1;                 % size of timeslice in units:
para.tdvp.maxExpMDim = 2*10^3;          % maximum allowed a_n*a_(n+1)*d_k
para.tdvp.rescaling = 0;                % turn on/off rescaling in TDVP
para.rescaling = para.tdvp.rescaling;
bondDim = [1 results.D{end} 1];

% estimate the actual dimension in expm() from MPS dimensions.
if para.useVmat
    para.tdvp.currentExpMDim = results.d_opt{end}.*bondDim(1:end-1).*bondDim(2:end);
else
    para.tdvp.currentExpMDim = results.dk{end}.*bondDim(1:end-1).*bondDim(2:end);
end



%% 1. sweep: condition the ground state MPS and prepare hamiltonian operator terms H(n) and K(n)
    % is already done in op.opstorage and op.hlrstorage
    %
    % also prepare other stuff

%% 2. time sweep

for timeslice = 1:(para.tdvp.tmax/para.tdvp.deltaT)
    para.sweepto = 'r';
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

    %% calculate time-dependent expectation values
end


end