
function [A, E] = minimizeE_onesiteA(op, B, Alaststep,para,sitej)

% ******************** one-site optimization ***********************
tol=para.eigs_tol;
DAl = size(op.Hleft, 1);
DAr = size(op.Hright, 1);
d = size(B, 2);
[ll,rr,dd] = size(Alaststep);
if para.parity~='n'
    nzAlast=Alaststep(para.Anzi{sitej});
end

% calculation of Heff
M = size(op.h2j, 1);

op.h1j = contracttensors(op.h1j,2,2,B,2,1);
op.h1j = contracttensors(conj(B),2,1,op.h1j,2,1);

for i=1:M
    op.h2j{i,1} = contracttensors(op.h2j{i,1},2,2,B,2,1);
    op.h2j{i,1} = contracttensors(conj(B),2,1,op.h2j{i,1},2,1);

    op.h2j{i,2} = contracttensors(op.h2j{i,2},2,2,B,2,1);
    op.h2j{i,2} = contracttensors(conj(B),2,1,op.h2j{i,2},2,1);
end


% projection on orthogonal subspace
%if ~isempty(P), Heff = P' * Heff * P; end

% optimization
opts.disp = 0;
opts.tol = tol;
%opts.p=min(24,DAl*DAr*d);
opts.issym =1;
sigma='sa';
if para.complex==1
    opts.isreal = false;
    sigma='sr';
end
assert(DAl*DAr*d == ll*rr*dd)

if para.parity=='n'
    opts.v0=reshape(Alaststep,ll*rr*d,1);
    [Avec, E]=eigs(@(x) HmultA(x, op, DAl, DAr, d, M,para.parity),DAl*DAr*d,1,sigma,opts);
    A = reshape(Avec, [DAl, DAr, d]);
else
    opts.v0=nzAlast;
    [Avec, E]=eigs(@(x) HmultA(x, op, DAl, DAr, d, M,para.parity,para.Anzi{sitej}),DAl*DAr*d/2,1,sigma,opts);
    A = zeros(DAl, DAr, d);
    A(para.Anzi{sitej})=Avec;
end

end
