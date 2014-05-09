function [c]=correlator_allsites(c_op,mps,Vmat)
% Calculate the expectation value of the multi-site correlator "c_op" for all bath sites. Only for successive operators. Handles any number M!
% Can be used instead of expectation_allsites if M = 1;
% c(j)=<\psi|c_op{j}|\psi>
% e.g.: c_op{j} = {bp{j} , bp{j+1}} == bp{j}*bp{j+1}*....*bp{j+m-1}
%	Modified:
%		FS 10/03/2014

% N=length(c_op);		%old definition, DELETE
[M, N] = size(c_op);											% N != Length, M = number of
assert(N==length(mps) && N==length(Vmat));

c=zeros(1,N);
cdset=cell(1,N);					% Used for complete representation including unity operators for one c_op{j};

for j=1:N
        cdset{1,j}=eye(size(Vmat{j},1));
end

for j=1:N
% 	for i = 1:M                                         % can be deleted
%         if j+i-1 <= N                                               % don't do for last site!
%             temp{i} = cdset{1,j+i-1};								% store the unity ops. TODO: might not be fast enough!
%             cdset{1,j+i-1}=c_op{i,j};								% produce: 1 .... 1 n 1 .... 1
%         end
% 	end
%
% This should be a bit faster than for-loop
    nterms = min(N-j+1,M);                                          % caps last terms for the end of the chain
    temp = cdset(1,j:min(j+M-1,N));                                 % save unity ops
    cdset(1,j:min(j+M-1,N)) = c_op(1:nterms,j)';                    % use c_ops
    c(j) = expectationvalue(cdset,mps,Vmat,mps,Vmat);               % calculate correlator
    cdset(1,j:min(j+M-1,N)) = temp;                                 % load back unity ops
% 	for i = 1:M                                         % can be deleted
%         if j+i-1 <= N
%     		cdset{1,j+i-1} = temp{i};								% write unity back.
%         end
% 	end
end