function [mps, Vmat, para, results, Hn] = tdvp_1site_evolveHn(mps,Vmat,para,results,op,sitej)
%% Evolves one site following Haegeman 2014
%   - Only contains A = exp(-i H(n) dt/2) A;
%   - Same procedure for l->r or l<-r
%   - Splits time-evolution into A and V part if para.useVmat == 1
%
% Created by Florian Schroeder @ Cambridge 20/10/2014
%

[BondDimLeft, BondDimRight, OBBDim]  = size(mps{sitej});
dk = results.dk{end}(sitej);
if ~para.useVmat
    assert(OBBDim == dk) ;
end
M = size(op.h2j,1);

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
    if size(HAA,1) <= para.tdvp.maxExpMDim
        Vmat_focused = expm(- 1i .* HAA .*para.tdvp.deltaT./2) * reshape(Vmat_focused,[dk*OBBDim,1]);
    else
        % Do approximation of exp(A)*v
        Vmat_focused = expv(- 1i*para.tdvp.deltaT./2,HAA,...
            reshape(Vmat_focused,[dk*OBBDim,1]),...
            para.tdvp.expvTol, para.tdvp.expvM);
    end
    Vmat_focused = reshape(Vmat_focused,[dk,OBBDim]);
    clear('HAA');

    %% normalise Vmat and take focus to A
    [Vmat{sitej}, V, results] = prepare_onesiteVmat(Vmat_focused,para,results,sitej);  % TODO: enable
%     [Vmat_focused, V, results] = prepare_onesiteVmat(Vmat_focused,para,results,sitej);  % TODO: disable
    mps{sitej} = contracttensors(Amat, 3, 3, V, 2, 2);     % TODO: enable later
%     Amat = contracttensors(Amat, 3, 3, V, 2, 2);            % TODO: disable
    clear('Amat','Vmat_focused','V');

end

%% Now: construct H(n)
% according to Haegeman 2014
% 1. Bring Hamiltonian terms into OBB if needed.

if para.useVmat     % contract H-terms to OBB; also ok for spinsites! since Vmat = eye
    % h1term to OBB, h1j can be rescaled
    % h1j_(n~',n~) = V*_(n',n~') [h1j_(n',n) V_(n,n~)]_(n',n~)
    h1j = contracttensors(op.h1j,2,2,Vmat{sitej},2,1);            % h1j_(n',n~)  = h1j_(n',n) V_(n,n~)
    h1j = contracttensors(conj(Vmat{sitej}),2,1,h1j,2,1);       % h1j_(n~',n~) = V*_(n',n~') h1j_(n',n~)

    % h2term to OBB, h2j can be rescaled
    % h2j_(n~',n~) = V*_(n',n~') [h2j_(n',n) V_(n,n~)]_(n',n~)
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
    t = para.tdvp.deltaT./2;
else
    t = para.tdvp.deltaT;
end

if size(Hn,1) <= para.tdvp.maxExpMDim
    mpsNew = expm(- 1i .* Hn .*t) * reshape(mps{sitej},[BondDimLeft*BondDimRight*OBBDim,1]);
else
    mpsNew = expv(- 1i*t, Hn,...
        reshape(mps{sitej},[BondDimLeft*BondDimRight*OBBDim,1]),...
        para.tdvp.expvTol, para.tdvp.expvM);
end

mps{sitej} = reshape(mpsNew,[BondDimLeft,BondDimRight,OBBDim]);
clear('mpsNew');
% now: A and V are time-evolved, A is focused
% if sitej = L, then start lr sweep with decomposition of mps


end