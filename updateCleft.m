function [Cleft] = updateCleft(Cleft, B, BUb, X, A, AUb)
    if isempty(X), X = speye(size(BUb, 1)); end
    if isempty(Cleft), Cleft=eye(size(B,1)); end

    newX = contracttensors(X,2,1,conj(BUb),2,1);
    newX = contracttensors(newX,2,1,AUb,2,1);


    Cleft = contracttensors(A, 3, 1, Cleft, 2, 2);
    Cleft = contracttensors(newX, 2, 2, Cleft, 3, 2);
    Cleft = contracttensors(conj(B), 3, [1, 3], Cleft, 3, [3, 1]);

