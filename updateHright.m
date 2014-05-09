function [Hright] = updateHright(Hright, h1j, Opstorage2, B, BUb, h2j, A, AUb, M)
% called as op.Hlrstorage{j} = updateHright(op.Hlrstorage{j + 1},  op.h1term{j}, op.Opstorage(:,2,j+1), mps{j}, Vmat{j},   op.h2term(:,1,j),  mps{j},  Vmat{j},M);	for middle sites	in initstorage.m

% As called from initstorage.m:
% Hright = op.Hlrstorage{j + 1}; == 0 for j=L
% h1j = op.h1term{j};	current on-site energy
% Opstorage2{}  = op.Opstorage(:,2,j+1); per site and  contains M matrices for each term in sum
% B = mps{j};
% BUb = Vmat{j};
% h2j = op.h2term(:,1,j);
% A = mps{j};
% AUb = Vmat{j}

% Commented by Florian Schroeder 29/01/2014

        if isempty(Opstorage2{1}) 				% initialize matrices if empty
            for m=1:M
                Opstorage2{m} = zeros(size(B,2));	% zeros(D)? (also for first site); zeors(1) for last
            end
        end

        if isempty(h2j{1})					% initialize matrices if empty
            for m=1:M
                h2j{m} = zeros(size(BUb,1));		% zeros(para.dk(j))
            end
        end


        Hright=updateCright(Hright, B, BUb, [], A, AUb);				% transform old Hright from effective basis j+1 into j
        if ~isempty(h1j)
            Hright = Hright + updateCright([],B,BUb,h1j,A,AUb);		% transform H1{j] into effective right basis of j-1
        end

        for m=1:M
            Hright = Hright+updateCright(Opstorage2{m},B,BUb,h2j{m},A,AUb);	% transform the interaction terms O_j O_{j+1} in the eff. right basis r_{j-1}
        end
end