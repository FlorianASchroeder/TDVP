function [Hleft] = updateHleft(Hleft, h1j, Opstorage1, B, BUb, h2j, A, AUb, M, varargin)
%% function [Hleft] = updateHleft(Hleft, h1jOBB, Opstorage1, B, BUb, h2jOBB, A, AUb, M, varargin)
%
%  Calculates the left effective Hamiltonian before moving to the next site
%  j+1.
%
% Modified
%	FS 17/07/2015: - added support for multi-chain on OBB level

	multiChain = 0;
	if nargin > 9
		para = varargin{1};
		if para.nChains > 1
			multiChain = 1;						% multiChain at OBB level
		end
	end

	if isempty(Opstorage1{1})
		for m=1:M
			Opstorage1 = zeros(size(B,1));
		end
	end

	if multiChain == 0 && isempty(h2j{1})
		for m=1:M
			h2j{m} = zeros(size(BUb,1));
		end
	end

	Hleft = updateCleft(Hleft, B, BUb, [], A, AUb);

	if ~isempty(h1j)
		Hleft = Hleft + updateCleft([],B,BUb,h1j,A,AUb);
	end

	for m=1:M
		Hleft = Hleft + updateCleft(Opstorage1{m},B,BUb,h2j{m},A,AUb);
	end


end