    function [Cright] = updateCright(Cright, B, BUb, X, A, AUb)

        if isempty(X), X = speye(size(BUb, 1)); end
        if isempty(Cright), Cright = eye (size(B,2)); end
        if isempty(BUb) && isempty(AUb)
            newX=X;
        else
        newX = contracttensors(X,2,1,conj(BUb),2,1);
        newX = contracttensors(newX,2,1,AUb,2,1);
        end

        Cright = contracttensors(A, 3, 2, Cright, 2, 2);
        Cright = contracttensors(newX, 2, 2, Cright, 3, 2);
        Cright = contracttensors(conj(B), 3, [2, 3], Cright, 3, [3, 1]);