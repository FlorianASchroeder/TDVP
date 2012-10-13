function [e] = expectationvalue(hset,mps1,V1,mps0,V0)
%Calculate the overlape <mps1|hset|mps0>
%If hset is empty then, it's eqvilent to calculate the overlap.

[M, N] = size(hset);

% expectation value
e = 0;
for m = 1:M
    em = 1;
    for j = N:-1:1
        em = updateCright(em, mps0{j},V0{j}, hset{m,j}, mps1{j},V1{j});
    end
    e = e + em;
end

% % norm
% n0=1;
% n1 = 1;
%
% for j = N:-1:1
%     d = size(mps0{j}, 3);
%     X = eye(d);
%    n0 = updateCright(n0, mps0{j}, [],X, mps0{j},[]);
% end
%
% for j = N:-1:1
%     d = size(mps1{j}, 3);
%     X = eye(d);
%    n1 = updateCright(n1, mps1{j}, [],X, mps1{j},[]);
% end
%
% e = e/sqrt(n1)/sqrt(n0);


end