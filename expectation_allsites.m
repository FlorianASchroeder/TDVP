function [n]=expectation_allsites(n_op,mps,Vmat)
%Calculate the expectation value of the local operater "n_op" for all sites.
% n(j)=<\psi|n_op{j}|\psi>

N=length(n_op);
assert(N==length(mps) && N==length(Vmat));

n=zeros(1,N);
ndset=cell(1,N);

for j=1:N
        ndset{1,j}=eye(size(Vmat{j},1));
end

for j=1:N
    temp=ndset{1,j};
    ndset{1,j}=n_op{j};
    n(j)=expectationvalue(ndset,mps,Vmat,mps,Vmat);
    ndset{1,j}=temp;
end