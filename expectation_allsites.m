function [n]=expectation_allsites(n_op,mps,Vmat)
% Calculate the expectation value of the local operater "n_op" for all sites.
% performs 2*3*N contractions
% assumes precise right-normalization
% n(j)=<\psi|n_op{j}|\psi>

N=length(n_op);
assert(N==length(mps) && N==length(Vmat));			% would also work for N < dim(mps)

n=zeros(1,N);

Cl = [];						% contains left part in effective j-basis
for j = 1:N
	% calculate site operator
	% j == 1 -> Cl = [];
	n(1,j) = trace(updateCleft(Cl,mps{j},Vmat{j},n_op{1,j},mps{j},Vmat{j}));	% this is the result!

	% take Cleft to next site
	Cl = updateCleft(Cl,mps{j},Vmat{j},[],mps{j},Vmat{j});
end

% old routine much slower! 3*N*N contractions
% ndset=cell(1,N);
%
% for j=1:N
%         ndset{1,j}=eye(size(Vmat{j},1));
% end
%
% for j=1:N
%     temp=ndset{1,j};
%     ndset{1,j}=n_op{j};											% produce: 1 .... 1 n 1 .... 1
%     n(j)=expectationvalue(ndset,mps,Vmat,mps,Vmat);
%     ndset{1,j}=temp;
% end