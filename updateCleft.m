function [Cleft] = updateCleft(Cleft, B, BUb, X, A, AUb)
% Lookup: Schollwoeck 2011 4.2.1: Contracts Left contraction residue of chain with current site A-matrices and Operators.
% 1. transforms Cleft from effective left basis of site k into effective left basis of k+1 if:
%		X = []	; B = A = mps{k}	; BUb = AUb = Vmat{k}
% 2. transforms operator X of site k into effective left basis of site k+1 if:
%		Cleft = []	; B = A = mps{k}	; BUb = AUb = Vmat{k}
% 3. Does 1 and 2 altogether if everything is given.

% Commented by Florian Schroeder 29/01/2014

	if isempty(X), X = speye(size(BUb, 1)); end
    if isempty(Cleft), Cleft=eye(size(B,1)); end

	% transform X into OBB using Vmat
    newX = contracttensors(X,2,1,conj(BUb),2,1);
    newX = contracttensors(newX,2,1,AUb,2,1);

    % do contraction:  C_fb = B*_dfe  Xnew_ec  A_abc  C_da	where 3rd indices are running over n_k
    Cleft = contracttensors(A, 3, 1, Cleft, 2, 2);                      % Cleft_bcd = A_abc C_da
    Cleft = contracttensors(newX, 2, 2, Cleft, 3, 2);                   % Cleft_ebd = Xnew_ec Cleft_bcd
    Cleft = contracttensors(conj(B), 3, [1, 3], Cleft, 3, [3, 1]);      % Cleft_fb = B*_dfe Cleft_ebd

