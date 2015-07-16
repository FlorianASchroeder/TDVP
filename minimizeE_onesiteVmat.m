% ******************** one-site optimization ***********************
function [B, E] = minimizeE_onesiteVmat(op, A, Blaststep,para)
% A is the local MPS A-matrix
% B is the optimized local optimal phonon basis transform matrix V
% Blaststep is the to be optimized Vmat matrix.

% Commented by Florian Schroeder 31/01/2014

tol=para.eigs_tol;
M = size(op.h2j, 1);
[Dblast,dlast]=size(Blaststep);							% dk x d_opt

op.HlOPB = contracttensors(op.Hleft,2,2,A,3,1);						% left non-interacting
op.HlOPB = contracttensors(conj(A),3,[1,2],op.HlOPB,3,[1,2]);

op.HrOPB = contracttensors(A,3,2,op.Hright,2,2);					% right non-interacting
op.HrOPB = contracttensors(conj(A),3,[1,2],op.HrOPB,3,[1,3]);

op.OpleftOPB= cell(M,1);
op.OprightOPB= cell(M,1);

for m=1:M
	op.OpleftOPB{m}= contracttensors(op.Opleft{m}, 2,2, A,3,1);		% left interacting
	op.OpleftOPB{m}= contracttensors(conj(A),3,[1,2],op.OpleftOPB{m},3,[1,2]);

	op.OprightOPB{m} = contracttensors(A,3,2,op.Opright{m},2,2);	% right interacting
	op.OprightOPB{m} = contracttensors(conj(A),3,[1,2],op.OprightOPB{m},3,[1,3]);
end


d = size(op.HlOPB, 1); 							% = d_opt
if iscell(op.h1j)
	Db = prod(cell2mat(cellfun(@(x) size(x,1),op.h1j, 'UniformOutput',false)));		% = dk for multi-chain Hamiltonians
else
	Db = size(op.h2j{1,1}, 1);                      % = dk
end

% projection on orthogonal subspace
%if ~isempty(P), Heff = P' * Heff * P; end

% optimization
opts.disp = 0;
opts.tol = tol;
%opts.p=12;
opts.issym =1;
sigma='sa';
if para.complex==1
    opts.isreal = false;
    sigma='sr';
end
assert(Db*d == Dblast*dlast);

if para.parity~='n'
    opts.v0=vertcat(reshape(Blaststep(1:Db/2,1:d/2),Db*d/4,1),reshape(Blaststep(Db/2+1:end,d/2+1:end),Db*d/4,1));
    [Bvec, E]=eigs(@(x) HmultVmat(x, op, Db, d, M, para.parity), Db*d/2,1,sigma,opts);
    B = zeros(Db, d);
    B(1:Db/2,1:d/2)=reshape(Bvec(1:Db*d/4),Db/2,d/2);
    B(Db/2+1:end,d/2+1:end)=reshape(Bvec(Db*d/4+1:end),Db/2,d/2);
else
    opts.v0=reshape(Blaststep,Db*d,1);
    [Bvec, E]=eigs(@(x) HmultVmat(x, op, Db,d, M,para.parity), Db*d,1,sigma,opts);
    B = reshape(Bvec, [Db, d]);
end
fprintf('\n Vmat E: %g',E);


end