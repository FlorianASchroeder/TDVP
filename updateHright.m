function [Hright] = updateHright(Hright, h1j, Opstorage2, B, BUb, h2j, A, AUb, M)

        if isempty(Opstorage2{1})
            for m=1:M
                Opstorage2{m} = zeros(size(B,2));
            end
        end

        if isempty(h2j{1})
            for m=1:M
                h2j{m} = zeros(size(BUb,1));
            end
        end


        Hright=updateCright(Hright, B, BUb, [], A, AUb);
        if ~isempty(h1j)
            Hright=Hright+updateCright([],B,BUb,h1j,A,AUb);
        end

        for m=1:M
            Hright = Hright+updateCright(Opstorage2{m},B,BUb,h2j{m},A,AUb);
        end
end