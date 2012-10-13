function [op] = initstorage(mps, Vmat,op,para)
% op.Hlrstorage stores the sub hamiltonian of the left and right block.
% op.Opstorage storse the operaters in the left and right bases.

Lam=para.Lambda;
L=length(mps);
M = size(op.h2term,1);
op.Hlrstorage = cell(1, L + 1);
op.Opstorage = cell(M, 2, L+1);

op.Hlrstorage{1, 1} = 0;
op.Hlrstorage{1, L + 1} = 0;

op.Hlrstorage{L} = updateHright(op.Hlrstorage{L + 1}, op.h1term{L}, op.Opstorage(:,2,L+1), mps{L}, Vmat{L}, op.h2term(:,1,L), mps{L}, Vmat{L},M);

for m=1:M
    op.Opstorage{m,2,L}=updateCright([],mps{L},Vmat{L},op.h2term{m,2,L},mps{L},Vmat{L});
    op.Opstorage{m,1,1}=0;
    op.Opstorage{m,2,L+1}=0;
end

for j = L-1:-1:2
    for m=1:M
        op.Opstorage{m,2,j}=updateCright([],mps{j},Vmat{j},op.h2term{m,2,j},mps{j},Vmat{j});
    end

    h2j=op.h2term(:,1,j);
    op.Hlrstorage{j} = updateHright(op.Hlrstorage{j + 1}, op.h1term{j}, op.Opstorage(:,2,j+1), mps{j}, Vmat{j},h2j, mps{j}, Vmat{j},M);
end
