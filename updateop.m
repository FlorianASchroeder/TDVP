function op=updateop(op,mps,Vmat,sitej,para)
%Update op when finished calculating at site sitej
%
% Changed:
%   FS 20/10/2014: - using switch instead of if-statement
%                  - removed op.Hright = op.Hlrstorage{sitej + 1} statement
%                    in case of scaling since this would always be the case
%                    then.
%	FS 16/07/2015: - added multi-chain support at OBB level!
%	FS 05/08/2015: - using h1jOBB, h2jOBB -> reverse previous changes!

M = para.M;
switch para.sweepto
    case 'r'
        % still use scaled op.Hleft to accomodate energy minimisation in
		if isempty(op.h1jOBB) || isempty(op.h1j)
			if ~iscell(Vmat{sitej})                                     % old VMPS and OBB-MultiChain code
				[op] = H_Eff([]  , Vmat{sitej}, 'A' , op, para);		% bring into OBB -> wraps all Multi-Chain magic
			else
				%% New V-Tensor MultiChain code
				[op] = H_Eff([]  , Vmat{sitej}, 'MC-OBB', op, para);
				[op] = H_Eff([]  , Vmat{sitej}, 'MC-A'  , op, para);
			end
		end
        % sweep to right already computed OBB
        op.Hlrstorage{sitej + 1}      = updateHleft(op.Hleft, op.h1jOBB, op.Opleft(:), mps{sitej}, Vmat{sitej}, op.h2jOBB(:,2), mps{sitej}, Vmat{sitej}, M, para);
		for m = 1:M
			if ~isempty(op.h2jOBB{m,1})
				op.Opstorage{m,1,sitej+1} = updateCleft([],mps{sitej},[],op.h2jOBB{m,1},mps{sitej},[]);			% replaced h2term by h2jOBB. is this right??
			else
				op.Opstorage{m,1,sitej+1} = [];
			end
		end
    case 'l'
		% if adjustdopt -> empty h1jOBB, h2jOBB -> recompute from h1term, h2term. Else, use them!
        % collect all j-1 < parts of Hamiltonian which are not interacting with j-1, and transforms into eff basis r_{j-1}
		if isempty(op.h1jOBB) || isempty(op.h1j)
			if ~iscell(Vmat{sitej})                                     % old VMPS and OBB-MultiChain code
				[op] = H_Eff([]  , Vmat{sitej}, 'A' , op, para);		% bring into OBB -> wraps all Multi-Chain magic
			else
				%% New V-Tensor MultiChain code
				[op] = H_Eff([]  , Vmat{sitej}, 'MC-OBB', op, para);
				[op] = H_Eff([]  , Vmat{sitej}, 'MC-A'  , op, para);
			end
		end
		op.Hlrstorage{sitej}		  = updateHright(op.Hlrstorage{sitej + 1}, op.h1jOBB, op.Opstorage(:,2,sitej+1),mps{sitej},Vmat{sitej}, op.h2jOBB(:,1), mps{sitej},Vmat{sitej}, M, para);
		for m = 1:M
			% 2nd term of interaction term into ef. basis of r_j-1; (1st is already in site basis.)
			if ~isempty(op.h2jOBB{m,2})
				op.Opstorage{m, 2, sitej} = updateCright([], mps{sitej},[], op.h2jOBB{m,2}, mps{sitej},[]);
			else
				op.Opstorage{m, 2, sitej} = [];
			end
		end
end

end