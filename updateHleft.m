function [Hleft] = updateHleft(Hleft, h1j, Opstorage1, B, BUb, h2j, A, AUb, M)

		if isempty(Opstorage1{1})
			for m=1:M
				Opstorage1 = zeros(size(B,1));
			end
		end

		if isempty(h2j{1})
			for m=1:M
				h2j{m} = zeros(size(BUb,1));
            end
		end

		Hleft=updateCleft(Hleft, B, BUb, [], A, AUb);

        if ~isempty(h1j)
			Hleft=Hleft+updateCleft([],B,BUb,h1j,A,AUb);
		end

		for m=1:M
            Hleft = Hleft+updateCleft(Opstorage1{m},B,BUb,h2j{m},A,AUb);
		end
end