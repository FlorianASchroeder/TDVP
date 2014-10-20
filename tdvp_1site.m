function tdvp_1site(mps,Vmat,para,results,op)
%% VMPS Time-dependent variational principle
%   implementation following Haegeman et al. 2014, arXiv:1408.5056
%
%   (only use without OBB (Vmat) at the moment!)
%
%   Created by Florian Schroeder 07/10/2014
%
%

%% 0. Parameter settings:

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

for timeslice = 1:tmax
    para.sweepto = 'r';
    for sitej = 1:para.L
        %% 3. Time evolve single Site!
        fprintf('%g-', sitej);
        op = gen_sitej_op(op,para,sitej,results.leftge);                     % take Site h1 & h2 Operators apply rescaling to Hleft, Hright, Opleft ...???

        BondDimLeft  = size(op.Hlrstorage{sitej},1);
        BondDimRight = results.D{end}(sitej);   % or size(op.Hlrstorage{sitej+1},1), if results.D is not properly written. But should work always since any last change to D is also last entry in that cell array!
        if para.useVmat
            OBBDim      = results.d_opt{end}(sitej);
        else
            OBBDim      = results.dk{end}(sitej);
        end
        dk           = results.dk{end}(sitej);

        %% If using Vmat, evolve it first, only BOSON!
        if para.useVmat == 1 && prod(sitej ~= para.spinposition)                % if bosonic site only!
            [Amat,V] = prepare_onesiteAmat(mps{sitej},para,sitej);              % left-normalize A, SVD in n.
            Vmat_focused = contracttensors(Vmat{sitej}, 2, 2, V, 2, 2);         % set focus on Vmat: V_(n,n~)
            clear('V');
            % Amat = MPS{sitej} left normalised;

            % now: Construct
            % HAA_(n',n~',n,n~) = H(n)_(l',r',n',l,r,n)*A*_(l',r',n~')*A_(l,r,n~)
            % for: V(t+dt) = exp(-i HAA dt)_(n',n~',n,n~) * V(t)_(n,n~)

            HAA = 0;    % kron( zeros(OBBDim), zeros(dk))

            %% non-interacting Hamiltonian terms: Hleft + Hmid + Hright
            % contracted with MPS
            %
            % Hleft_(n~',n~) = A*_(l',r,n~') [Hl_(l',l) * A_(l,r,n~)]_(l',r,n~)
            Hleft = contracttensors(op.Hleft,2,2,Amat,3,1);                     % Hleft_(l',r,n~) = Hl_(l',l) * A_(l,r,n~)
            Hleft = contracttensors(conj(Amat),3,[1 2],Hleft,3,[1 2]);          % Hleft_(n~',n~)  = A*_(l',r,n~') * Hleft_(l',r,n~)
            HAA = HAA + kron(Hleft, eye(dk));
            clear('Hleft');

            % Hright_(n~',n~) = A*_(l,r',n~') [Hr_(r',r) * A_(l,r,n~)]_(r',l,n~)
            Hright = contracttensors(op.Hright,2,2, Amat,3,2);                  % Hright_(r',l,n~) = Hr_(r',r) * A_(l,r,n~)
            Hright = contracttensors(conj(Amat),3,[1 2], Hright,3,[2 1]);       % Hright_(n~',n~)  = A*_(l',r',n~') * Hright_(r',l,n~)
            HAA = HAA + kron(Hright, eye(dk));
            clear('Hright');

            % Hmid_(n~',n~) = A*_(l,r,n~') * A_(l,r,n~)
            Hmid = contracttensors(conj(Amat),3,[1 2], Amat,3,[1 2]);           % should be eye(OBBDim) TODO: could be replaced.
            HAA = HAA + kron(Hmid, op.h1j);
            clear('Hmid');

            %% Interacting Hamiltonian terms: \sum_i^M op.Opleft
            for m = 1:M
                % Opleft_(n~',n~) = A*_(l',r,n~') [Opleft_(l',l) * A_(l,r,n~)]_(l',r,n~)
                Opleft = contracttensors(op.Opleft{m},2,2,Amat,3,1);            % Opleft_(l',r,n~) = Opleft_(l',l) * A_(l,r,n~)
                Opleft = contracttensors(conj(Amat),3,[1 2],Opleft,3,[1 2]);    % Opleft_(n~',n~)  = A*_(l',r,n~') * Opleft_(l',r,n~)
                HAA = HAA + kron(Opleft, op.h2j{m,2});
                clear('Opleft');

                % Opright_(n~',n~) = A*_(l,r',n~') [Opright_(r',r) * A_(l,r,n~)]_(r',l,n~)
                Opright = contracttensors(op.Opright{m},2,2, Amat,3,2);         % Opright_(r',l,n~) = Opright_(r',r) * A_(l,r,n~)
                Opright = contracttensors(conj(Amat),3,[1 2], Opright,3,[2 1]); % Opright_(n~',n~)  = A*_(l',r',n~') * Opright_(r',l,n~)
                HAA = HAA + kron(Opright, op.h2j{m,1});
                clear('Opright');
            end

            %% Take matrix exponential
            % V(t+dt) = exp(-i HAA dt)_(n',n~',n,n~) * V(t)_(n,n~)
            HAAexp = expm(- 1i .* HAA .*para.tdvp.deltaT./2);
            Vmat_focused = HAAexp * reshape(Vmat_focused,[dk*OBBDim,1]);
            Vmat_focused = reshape(Vmat_focused,[dk,OBBDim]);
            clear('HAAexp','HAA');

            %% normalise Vmat and take focus to A
            %[Vmat{sitej}, V, results] = prepare_onesiteVmat(Vmat_focused,para,results,sitej);  % TODO: enable
            [Vmat_focused, V, results] = prepare_onesiteVmat(Vmat_focused,para,results,sitej);  % TODO: disable
            %mps{sitej} = contracttensors(Amat, 3, 3, V, 2, 2);     % TODO: enable later
            Amat = contracttensors(Amat, 3, 3, V, 2, 2);            % TODO: disable

        end

        %% Now: construct H(n)
        % according to Haegeman 2014
        % 1. Bring Hamiltonian terms into OBB if needed.

        if para.useVmat     % contract H-terms to OBB
            % h1term to OBB, h1j can be rescaled
            % h1j_(n~',n~) = V*_(n',n~') [h1j_(n',n) V_(n,n~)]_(n',n~)
            h1j = contracttensors(op.h1j,2,2,Vmat{sitej},2,1);            % h1j_(n',n~)  = h1j_(n',n) V_(n,n~)
            h1j = contracttensors(conj(Vmat{sitej}),2,1,h1j,2,1);       % h1j_(n~',n~) = V*_(n',n~') h1j_(n',n~)

            % h2term to OBB, h2j can be rescaled
            % h2j_(n~',n~) = V*_(n',n~') [h2j_(n',n) V_(n,n~)]_(n',n~)
            M = size(op.h2j, 1);
            h2j = cell(M,2);
            for i=1:M
                h2j{i,1} = contracttensors(op.h2term{i,1,sitej},2,2,Vmat{sitej},2,1);
                h2j{i,1} = contracttensors(conj(Vmat{sitej}),2,1,h2j{i,1},2,1);

                h2j{i,2} = contracttensors(op.h2term{i,2,sitej},2,2,Vmat{sitej},2,1);
                h2j{i,2} = contracttensors(conj(Vmat{sitej}),2,1,h2j{i,2},2,1);
            end
        else                % no OBB, then OBBDim = dk
            h1j = op.h1j;
            h2j = op.h2j;
        end

        Hn=0;       % Hn = kron(eye(OBBDim),kron(eye(BondDimRight),eye(BondDimLeft)))
        % all terms:
        Hn = Hn + kron(eye(OBBDim),kron(eye(BondDimRight),op.Hleft));
        Hn = Hn + kron(eye(OBBDim),kron(op.Hright,eye(BondDimLeft)));
        Hn = Hn + kron(h1j,kron(eye(BondDimRight),eye(BondDimLeft)));
        for m=1:M
            Hn = Hn + kron(h2j{m,2},kron(eye(BondDimRight),op.Opleft{m}));
            Hn = Hn + kron(h2j{m,1},kron(op.Opright{m},eye(BondDimLeft)));
        end

        %% Take and apply Matrix exponential
        % A(t+dt) = exp(-i Hn dt)_(l',r',n',l,r,n) * A(t)_(l,r,n)
        % Last site special, see Haegeman 2014
        % TODO: change expm() with threshold
        if sitej ~= para.L
            Hnexp = expm(- 1i .* Hn .*para.tdvp.deltaT./2);
        else
            Hnexp = expm(- 1i .* Hn .*para.tdvp.deltaT);
        end
        mpsNew = Hnexp * reshape(mps{sitej},[BondDimLeft*BondDimRight*OBBDim,1]);
        mps{sitej} = reshape(mpsNew,[BondDimLeft,BondDimRight,OBBDim]);
        clear('Hnexp','mpsNew');
        % now: A and V are time-evolved.
        % if sitej = L, then start lr sweep with decomposition of mps

        %% Take focus to next site
        if sitej ~= para.L
            %% Left-normalize A and get Center matrix C(n,t+dt)_(rl,r)
            [mps{sitej}, Cn, para,results] = prepare_onesite(mps{sitej}, para,sitej,results);

            %% prepare K(n) from H(n)
            % through contraction with A(t+dt)
            % K(n)_(rl',r',rl,r) = [A*_(l',rl',n') * H(n)_(l',r',n',l,r,n)]_(rl',r',l,r,n) * A_(l,rl,n)
            Hn = reshape(Hn, [BondDimLeft, BondDimRight, OBBDim, BondDimLeft, BondDimRight, OBBDim]);
            Kn = contracttensors(conj(mps{sitej}),3,[1 3], Hn,6,[1 3]);
            Kn = contracttensors(Kn,5,[3 5], mps{sitej},3,[1 3]);                   % K(n)_(rl',r',r,rl)
            Kn = permute(Kn,[1,2,4,3]);                                             % K(n)_(rl',r',rl,r)
            Kn = reshape(Kn,[BondDimRight^2, BondDimRight^2]);                      % K(n)_(rl'r',rlr)

            %% Take and apply Matrix exponential
            % C(n,t) = exp(+ i K(n) dt)_(rl'*r',rl*r) * C(n,t+dt)_(rl*r)
            Cn = reshape(Cn,[BondDimRight^2,1]);
            Cn = expm( 1i .* Kn .* para.tdvp.deltaT./2) *Cn;
            Cn = reshape(Cn, [BondDimRight, BondDimRight]);
            clear('Kn', 'Hn');

            %% Multiply Ac(n+1) = C(n) * Ar(n+1)
            % set focus on next site A
            % A_(rl,r,n) = C(n)_(rl,l) * A_(l,r,n)
            mps{sitej+1} = contracttensors(Cn,2,2, mps{sitej+1},3,1);
            clear('Cn');

            %% update Left Hamiltonian operators
            op = updateop(op,mps,Vmat,sitej,para);
        end
    end

    %% Log vNE etc ?
    % from prepare_onesite() and prepare_onesiteVmat():
    % results.Amat_vNE  array
    % results.Amat_sv   cell
    % results.Vmat_vNE
    % results.Vmat_sv

    %% turn around
    para.sweepto = 'l';

    for sitej = para.L-1:-1:1
        %% Right-normalize A(n+1) and get Center matrix C(n,t+dt)_(rl,r)
        [mps{sitej+1}, Cn, para,results] = prepare_onesite(mps{sitej+1},para,sitej,results);

        %% update right Hamiltonian operators
        op = updateop(op,mps,Vmat,sitej,para);
    end
end

%% test: Benchmark expm performance
finish = zeros(9,100);
for j=1:100
for i = 1:9
   a = randn(round(exp(i)));
   start = tic;
   expm(a);
   finish(i,j) = toc(start);
end
end


end