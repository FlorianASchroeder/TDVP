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
%       'spin':             cell, fields: sx, sy, sz
%       'occupation':       array
%       'shift':            array
%       'rdm':              matrix
%       'participation':    number
%       'tunnelenergy':     number
%
%
%   Created 03/06/2014 by Florian Schroeder @ University of Cambridge
%
%
switch type{1}
    case 'spin'
        % applicable for spin-boson model and for folded SBM2
        out = 0;

    case 'occupation'
        % applicable to all single stranded chains. Extensible to 2-chains
        out = 0;

    case 'shift'
        % applicable to all single stranded chains. Extensible to 2-chains
        out = 0;

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
end

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

reducedDensity = contracttensors(mps{k},3,[1,2],conj(mps{k}),3,[1,2]);  % contract rD_nm = Mk_abn Mk*_abm
reducedDensity = contracttensors(reducedDensity,2,2,conj(Vmat{k}),2,2);       % contract rD_nj = rD_nm Vmat*_jm
reducedDensity = contracttensors(Vmat{k},2,2,reducedDensity,2,1);       % contract rD_ij = Vmat_in rd_nj

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