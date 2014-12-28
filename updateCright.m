function [Cright] = updateCright(Cright, B, BUb, X, A, AUb)
% Computes contraction to put operator X in effective basis of B and A.
% BUb and AUb are the V matrices to express A in optimal boson basis.
% B and A are left and right MPS A-matrices
%
% 1. transforms Cright from effective right basis of site k into effective right basis of k-1 if:
%		X = []	; B = A = mps{k}	; BUb, AUb = Vmat{k}
% 2. transforms operator X of site k into effective right basis of site k-1 if:
%		Cleft = []	; B = A = mps{k}	; BUb = AUb = Vmat{k}
% 3. Does 1 and 2 altogether if everything is given.
%
% Called as updateCright([],mps{L},Vmat{L},op.h2term{m,2,L},mps{L},Vmat{L});				--> express in effective basis
% Called as updateCright(Hright, mps{j}, Vmat{j}, [], mps{j}, Vmat{j}); from updateHright.m		-->
%		Hright = op.Hlrstorage{j + 1} which == 0 for j=L
%
% if Cright input == 0 --> output == 0
% Commented by Florian Schroeder 29/01/2014
%
% Modified
%	- 21/12/14 FS: replaced OBB contraction by faster matrix product

        if isempty(X), X = speye(size(BUb, 1)); end 			% newX = eye = X, as Vmat unitary
        if isempty(Cright), Cright = eye (size(B,2)); end
        if isempty(BUb) && isempty(AUb)					% if no Vmat
            newX=X;
        else
			% do (Vmat^†) . X . Vmat.
			% express X in OBB by applying Vmat
			newX = BUb' * X * AUb;
% 			newX = contracttensors(X,2,1,conj(BUb),2,1);			% newX_nk = X_mn * conj(BUb)_mk = T(X)_nm * conj(BUb)_mk
% 			newX = contracttensors(newX,2,1,AUb,2,1);				% newX = newX_nk * AUb_nl = adj(BUb)_km * X_mn * AUb_nl;
		end

    % if Cright = eye: contract X (in OBB)  with A matrices to transform into effective basis representation.
	% do contraction:  C_fa = B*_fde  Xnew_ec  A_abc  C_db	where 3rd indices are running over n_k
        Cright = contracttensors(A, 3, 2, Cright, 2, 2);
        Cright = contracttensors(newX, 2, 2, Cright, 3, 2);
        Cright = contracttensors(conj(B), 3, [2, 3], Cright, 3, [3, 1]);	% Contracts 2-3 and 3-1. I think..