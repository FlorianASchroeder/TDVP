function [op] = initstorage(mps, Vmat,op,para)
% op.Hlrstorage stores the sub hamiltonian of the (left?) and right block.
% op.Opstorage stores the operators in the (left?) and right bases.
% sweeps from right to left.

% Commented by Florian Schroeder 29/01/2014

Lam=para.Lambda;			% unused! can be deleted
L=length(mps);
M = size(op.h2term,1); 		% number of terms in sum per site
op.Hlrstorage = cell(1, L + 1);
op.Opstorage = cell(M, 2, L+1);

% treat first and last terms seperately
op.Hlrstorage{1, 1} = 0;
op.Hlrstorage{1, L + 1} = 0;

op.Hlrstorage{L} = updateHright(op.Hlrstorage{L + 1}, op.h1term{L}, op.Opstorage(:,2,L+1), mps{L}, Vmat{L}, op.h2term(:,1,L), mps{L}, Vmat{L},M);
	% gives h1term{L} in effective basis of r_L-1 as:
	% op.Hlrstorage{l+1} = 0; op.h1term{L} = exists ; op.Opstorage(:,2,L+1) = [], op.h2term(:,1,L) = 0;

for m=1:M
    op.Opstorage{m,2,L}=updateCright([],mps{L},Vmat{L},op.h2term{m,2,L},mps{L},Vmat{L});
	% transforms interaction terms into r_L-1 ef. basis
    op.Opstorage{m,1,1}=0;
    op.Opstorage{m,2,L+1}=0;
end

% middle terms builds from r->l
for j = L-1:-1:2
    for m=1:M
        op.Opstorage{m,2,j}=updateCright([],mps{j},Vmat{j},op.h2term{m,2,j},mps{j},Vmat{j});
		% 2nd term of interaction term into ef. basis of r_j-1; (1st is already in site basis.)
    end

    h2j=op.h2term(:,1,j);
    op.Hlrstorage{j} = updateHright(op.Hlrstorage{j + 1}, op.h1term{j}, op.Opstorage(:,2,j+1), mps{j}, Vmat{j},h2j, mps{j}, Vmat{j},M);
	% collects all j-1 < parts of Hamiltonian which are not interacting with j-1, and transforms into eff basis r_{j-1}
end
