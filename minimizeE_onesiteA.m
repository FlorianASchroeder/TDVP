function [A, E] = minimizeE_onesiteA(op, B, Alaststep,para)
% B can be Vmat{sitej}; for spinsite: = eye()


% Commented by Florian Schroeder 31/01/2014
% Modified:
%	- FS 01/12/2014: replaced OBB transformation by matrix products.
%					 10x faster! A' = conjugate transpose!

% ******************** one-site optimization ***********************
tol=para.eigs_tol;
DAl = size(op.Hleft, 1);					% Bond dimension Dleft of A to left
DAr = size(op.Hright, 1);					% Bond dimension Dright of A to right
d = size(B, 2);								% = d_opt of Vmat
[ll,rr,dd] = size(Alaststep);				% = Dleft, Dright, d_opt	from A
if para.parity~='n'
    nzAlast=Alaststep(para.Anzi{para.sitej});
end

% calculation of Heff
M = size(op.h2j, 1);

% transform all bare H terms of sitej into OBB
if para.nChains == 1
	op.h1j = B' * (op.h1j * B);									% faster and more accurate
	% op.h1j = contracttensors(op.h1j,2,2,B,2,1);             % = h1_ab V_bd
	% op.h1j = contracttensors(conj(B),2,1,op.h1j,2,1);       % -> h1'_cd = V*_ac h1_ab V_bd
	for i=1:M
		op.h2j{i,1} = B' * (op.h2j{i,1} * B);									% faster and more accurate
	%     op.h2j{i,1} = contracttensors(op.h2j{i,1},2,2,B,2,1);
	%     op.h2j{i,1} = contracttensors(conj(B),2,1,op.h2j{i,1},2,1);
		op.h2j{i,2} = B' * (op.h2j{i,2} * B);									% faster and more accurate
	%     op.h2j{i,2} = contracttensors(op.h2j{i,2},2,2,B,2,1);
	%     op.h2j{i,2} = contracttensors(conj(B),2,1,op.h2j{i,2},2,1);
	end
else  % nChains > 1
	h1jnew = 0;
	for i = find(~cellfun('isempty',op.h1j'))
		H1 = cell(para.nChains,1);
		H1(i) = op.h1j(i);
		h1jnew = h1jnew + contractMultiChainOBB(B, H1, para);
	end
	op.h1j = h1jnew;

	h2jnew = cell(M,2);
	for i = 1:M
		h2jnew{i,1} = contractMultiChainOBB(B, op.h2j(i,1,:), para);
		h2jnew{i,2} = contractMultiChainOBB(B, op.h2j(i,2,:), para);
	end
	op.h2j = h2jnew;
end



% projection on orthogonal subspace
%if ~isempty(P), Heff = P' * Heff * P; end

%% optimization

% options for eigs()
opts.disp = 0;
opts.tol = tol;
%opts.p=min(24,DAl*DAr*d);
opts.issym =1;
sigma='sa';							% smallest algebraic eigenvalue
if para.complex==1
    opts.isreal = false;
    sigma='sr';						% smallest real part eigenvalue
end
assert(DAl*DAr*d == ll*rr*dd)

if para.parity=='n'
    opts.v0=reshape(Alaststep,ll*rr*d,1);						% starting guess, old MPS A-matrix
    [Avec, E]=eigs(@(x) HmultA(x, op, DAl, DAr, d, M,para.parity),DAl*DAr*d,1,sigma,opts);
    A = reshape(Avec, [DAl, DAr, d]);
else
    opts.v0=nzAlast;
    [Avec, E]=eigs(@(x) HmultA(x, op, DAl, DAr, d, M,para.parity,para.Anzi{para.sitej}),DAl*DAr*d/2,1,sigma,opts);
    A = zeros(DAl, DAr, d);
    A(para.Anzi{para.sitej})=Avec;
end
fprintf('\n Amat E: %g',E);

end
