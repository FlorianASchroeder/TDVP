function out = getObservable(type,mps,Vmat,para)
%% Calculates the following Observables:
%   SBM:    'spin', 'occupation', 'shift', 'rdm'
%   MLSBM:  'participation', 'tunnelenergy'
%
%   Use as: first argument is cell, with descriptor as type{1}
%               all other type{n} are options defined in each function
%           getObservable({'spin'},mps,Vmat,para)
%           getObservable({'rdm',2},mps,Vmat,para)
%           getObservable({'tunnelenergy',op},mps,Vmat,para)
%
%   The output 'out' can be of different form:
%       'spin':             struct, fields: sx, sy, sz
%       'occupation':       array
%       'shift':            array
%       'rdm':              matrix
%       'participation':    number
%       'tunnelenergy':     number
%
%
%   Created 03/06/2014 by Florian Schroeder @ University of Cambridge
%   TODO:   - Implement Boson Site Shift to export the routine from optimizesite.m. Only get a single shift value.
%
% Modified:
%	- 21/12/14 FS: replaced OBB contractions with faster matrix products.
switch type{1}
    case 'spin'
        % applicable for spin-boson model and for folded SBM2
        out = calSpin(mps,Vmat,para);

    case 'occupation'
        % applicable to all single stranded chains. Extensible to 2-chains
        out = calBosonOcc(mps,Vmat,para);

    case 'shift'
        % applicable to all single stranded chains. Extensible to 2-chains
        out = calBosonShift(mps,Vmat,para);

    case 'rdm'
        % applicable to all systems.
        if type{2} <= para.L
            out = calRDM(mps,Vmat,para,type{2});
        else
            out = [];
        end

    case 'participation'
        % applicable to all systems, only for first site.
        out = calParticipation(calRDM(mps,Vmat,para,1));

    case 'tunnelenergy'
        % useful for MLSBM, applicable to all
        % needs type{2} = op
        if length(type) == 2
            out = calTunnelingEnergy(mps,Vmat,para,type{2});
        else
            out = calTunnelingEnergy(mps,Vmat,para);
        end
    case 'coherence'
        % only for spin in SBM

end

end

function spin = calSpin(mps,Vmat,para)
% Calculate the spin expectation value
% Has to be modified if site 1 changes dimension!!
% Modified:
%   FS 24/05/2014:  excluded MLSpinBoson. Use with try-catch block
%
%   TODO: Needs extension to use information about dimension of spin site.
%

N=para.L;
assert(N==length(mps) && N==length(Vmat));
assert(~strcmp(para.model,'MLSpinBoson'),'not possible for MLSBM');

ndset=cell(1,N);

for j=1:N
    ndset{1,j}=eye(size(Vmat{j},1));
end
sx=ndset;sy=sx;sz=sy;

%debug:
% sx{1,1}

[sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
sx{para.spinposition}=sigmaX;
sy{para.spinposition}=sigmaY;
sz{para.spinposition}=sigmaZ;
if strcmp(para.model,'2SpinPhononModel')
    sx{para.spinposition}=kron(sigmaZ,eye(2));  %measures excitation of site1
    sz{para.spinposition}=kron(eye(2),sigmaZ);  %measures excitation of site2
    %try to find a good way to measure this!
    sy{para.spinposition}=kron(sigmaY,eye(2));  %measures only sy of site1
end
spin.sx=expectationvalue(sx,mps,Vmat,mps,Vmat);
spin.sy=expectationvalue(sy,mps,Vmat,mps,Vmat);
spin.sz=expectationvalue(sz,mps,Vmat,mps,Vmat);
spin.sx=real(spin.sx);
spin.sy=real(spin.sy);
spin.sz=real(spin.sz);

end

function nx = calBosonOcc(mps,Vmat,para)
% Calculate the boson occupation on x chain
% The operator on the spin site is set to zero
% modify this to get left and right chain occupation! perhaps in calbosonocc_SBM2.m
%
% Modified:
%	FS 23/01/2014:	- Introduced '~' to ignore unused returned values
%					- support for folded Chain models
%   FS 10/03/2014:  - Using correlator_allsites which is more general
%   FS 04/06/2014:  - updated for spinposition array.
%

n_op = cell(1,para.L);

for j=1:para.L
    if prod(j~=para.spinposition)
        if para.foldedChain == 0
            % 1-chain SBM
            [~,~,n] = bosonop(para.dk(j),para.shift(j),para.parity);
            n_op{j} = n;
        %Modification for 2chain model!! Not perfect or right yet!
        elseif para.foldedChain == 1
            % only kron(n,1) chain occupation calculated.
            if para.parity == 'n'
                [~,~,n] = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
                idm = eye(size(n));
                nr = kron(n,idm);
            else
                [bp,~,~] = bosonop(para.dk(j),para.shift(j),para.parity);   % Why without sqrt??
                [~,~,nr,~,~,~]=paritykron(bp,para.bosonparity);
            end
            n_op{j} = nr;
        else
            disp('Not Implemented: getObservable, calBosonOcc, foldedchain >1');
        end
    else
        n_op{j}=zeros(para.dk(j));      % don't measure spin.
    end
end

%
nx = correlator_allsites(n_op,mps,Vmat);

end

function bosonshift = calBosonShift(mps,Vmat,para)
% Calculate boson shift x, x^2, var(x)
% The operator on the spin site is set to zero
% Modified:
%	FS 22/01/2014:	- changed to using para.foldedChain.
%   FS 04/06/2014:  - updated for spinposition array
%
% expectation_allsites is old and should be replaced by correlator_allsites

x_opx=cell(1,para.L);
x2_opx=cell(1,para.L);
for j=1:para.L
    if prod(j~=para.spinposition)
        if para.foldedChain == 1
            % Constructs Supersite Operators
            % Only measures kron(x,1) chain part
            [bp,~,n] = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
            idm = eye(size(n));
            bpx = kron(bp,idm); bmx = bpx';         %nx = kron(n,idm); unused
            x_opx{j}  = sqrt(2)/2.*(bpx+bmx);
            x2_opx{j} = x_opx{j}*x_opx{j};
        elseif para.foldedChain == 0
            [bp,bm,~] = bosonop(para.dk(j),para.shift(j),para.parity);
            x_opx{j}  = sqrt(2)/2.*(bp+bm);
            x2_opx{j} = x_opx{j}*x_opx{j};
        end
    else
        x_opx{j}  = zeros(para.dk(j));
        x2_opx{j} = zeros(para.dk(j));
    end
end
bosonshift.x        = expectation_allsites(x_opx,mps,Vmat);
bosonshift.xsquare  = expectation_allsites(x2_opx,mps,Vmat);
bosonshift.xvariant = sqrt(bosonshift.xsquare-bosonshift.x.^2);
bosonshift.xerror   = mean(abs(para.shift-bosonshift.x));
end

function reducedDensity = calRDM(mps,Vmat,para,k)
% calculates the reduced density matrix of any single site, mps{k} for size(mps{k})=a_{k-1} x a_k x n_k
%
%   created 24/05/2014 by Florian Schroeder @ University of Cambridge
%
%
% copied from prepare.m:
% does l -> r sweep to create state in local picture of k

for i = 1:k-1
    if para.useVmat==1
        [Vmat{i},V] = prepare_onesiteVmat(Vmat{i},para);			% Vmat = U * S * V' ; Vmat := U; V:= S * V'
        mps{i} = contracttensors(mps{i},3,3,V,2,2);                 % = Ai_{l,r,n} * V'_{p,n}; This contraction is defined differently to the paper.
    end
    [mps{i}, U] = prepare_onesite(mps{i}, 'lr',para,i);             % SVD(Ai_(l,r,n)) = Ai_(l,m,n) * U_(m,r)
    mps{i+1} = contracttensors(U,2,2,mps{i+1},3,1);                 % U_(m,l) * A(i+1)_(l,r,n)
    para=gennonzeroindex(mps,Vmat,para,i);                          % only if parity not 'n'
    para=gennonzeroindex(mps,Vmat,para,i+1);                        % only if parity not 'n'
end

% now in form: Al{1}...Al{k-1} M{k} Ar{k+1}...Ar{L}
%   with Al = left-normalized, Ar: right-normalized.

reducedDensity = contracttensors(mps{k},3,[1,2],conj(mps{k}),3,[1,2]);		% contract rD_nm = Mk_abn Mk*_abm

% reducedDensity = Vmat{k} * (reducedDensity * Vmat{k}');
reducedDensity = contracttensors(reducedDensity,2,2,conj(Vmat{k}),2,2);     % contract rD_nj = rD_nm Vmat*_jm
reducedDensity = contracttensors(Vmat{k},2,2,reducedDensity,2,1);			% contract rD_ij = Vmat_in rd_nj

end

function participation = calParticipation(rdm)
% rdm: a reduced density matrix of a single site
%   e.g. rdm = calRDM(mps,Vmat,para,k)
%
% takes a rdm and calculates the participation ratio

participation = 1/sum(diag(rdm).^2);

end

function tunnelE = calTunnelingEnergy(mps,Vmat,para,op)
%% calculates the tunneling energy in MLSBM of the PPC system
%   <Psi|H0-diag(H0)|Psi>
%
%   could be applied to any system.
%   Interacting System Hamiltonian: HI


HI = cell(1);
if ~exist('op','var') || ~isfield(op,'h1term')
    disp('guessing H0 from para for MLSBM');
    op.h1term{1,1} = Hamiltonian_PPC(para);
end
HI{1} = op.h1term{1,1}-diag(diag(op.h1term{1,1}));
tunnelE = expectationvalue(HI,mps,Vmat,mps,Vmat);

end