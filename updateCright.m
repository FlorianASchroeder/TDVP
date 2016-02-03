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
	skipX = 0;
	if isempty(X) && isempty(BUb) && isempty(AUb)
		skipX = 1;
	elseif isempty(X) && ~isempty(BUb)
		X = speye(size(BUb, 1));					% newX = eye = X, as Vmat unitary
	end
	
	if isempty(BUb) && isempty(AUb)					% if no Vmat
		newX=X;
	else
		% do (Vmat^†) . X . Vmat.
		% express X in OBB by applying Vmat
		newX = (BUb' * X) * AUb;
	end

    % if Cright = eye: contract X (in OBB)  with A matrices to transform into effective basis representation.
	% do contraction:  C_fa = B*_fde  A_abc  Xnew_ec  C_db 	where 3rd indices are running over n_k
	if ~skipX
		A      = contracttensors(A, 3, 3, newX.', 2, 1);			% Anew_abe = A_abc  Xnew_ec
	end
	if ~isempty(Cright)
		A = contracttensors(A, 3, 2, Cright.', 2, 1);				% Anew_aed = A_abe C_db
		Cright = contracttensors(conj(B), 3, [2, 3], A, 3, [3, 2]);	% Contracts 2-3 and 3-2
	else
		%have Anew_ade since C_db = 1_db -> neglected
		Cright = contracttensors(conj(B), 3, [2, 3], A, 3, [2, 3]);	% Contracts 2-2 and 3-3
	end




