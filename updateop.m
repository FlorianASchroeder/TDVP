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
        % ground state search
		op.Hlrstorage{sitej + 1}      = updateHleft(op.Hleft, op.h1j, op.Opleft(:), mps{sitej}, Vmat{sitej}, op.h2j(:,2), mps{sitej}, Vmat{sitej}, M);
		for m = 1:M
			op.Opstorage{m,1,sitej+1} = updateCleft([],mps{sitej},Vmat{sitej},op.h2term{m,1,sitej},mps{sitej},Vmat{sitej});
		end
    case 'l'
		% if adjustdopt -> empty h1jOBB, h2jOBB -> recompute from h1term, h2term. Else, use them!
        % collect all j-1 < parts of Hamiltonian which are not interacting with j-1, and transforms into eff basis r_{j-1}
		op.Hlrstorage{sitej}		  = updateHright(op.Hlrstorage{sitej + 1}, op.h1term{sitej}, op.Opstorage(:,2,sitej+1),mps{sitej},Vmat{sitej}, op.h2term(:,1,sitej), mps{sitej},Vmat{sitej}, M);
		for m = 1:M
			% 2nd term of interaction term into ef. basis of r_j-1; (1st is already in site basis.)
			op.Opstorage{m, 2, sitej} = updateCright([], mps{sitej},Vmat{sitej}, op.h2term{m,2,sitej}, mps{sitej},Vmat{sitej});
		end
end
end