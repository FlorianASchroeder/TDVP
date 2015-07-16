function [Hleft] = updateHleft(Hleft, h1j, Opstorage1, B, BUb, h2j, A, AUb, M, varargin)
%% function [Hleft] = updateHleft(Hleft, h1j, Opstorage1, B, BUb, h2j, A, AUb, M, varargin)
%
%  Calculates the left effective Hamiltonian before moving to the next site
%  j+1.
%
% Modified
%	FS 17/07/2015: - added support for multi-chain on OBB level

	if nargin > 9
		para = varargin{1};
		multiChain = 1;						% multiChain at OBB level
	else
		multiChain = 0;
	end

	if isempty(Opstorage1{1})
		for m=1:M
			Opstorage1 = zeros(size(B,1));
		end
	end

	if multiChain == 0	&& any(cellfun('isempty',h2j))				% initialize matrices as 0 if empty (only use for single chains)
		for m = find(cellfun('isempty',h2j))
			h2j{m} = zeros(size(BUb,1));							% zeros(para.dk(j)) since empty h2j in updateCright -> eye -> wrong! need explicit 0!
		end															% multiChain == 1 does not use updateCright -> leave empty
	end

	Hleft=updateCleft(Hleft, B, BUb, [], A, AUb);

	if multiChain == 0 && ~isempty(h1j)
		Hleft=Hleft+updateCleft([],B,BUb,h1j,A,AUb);
	elseif multiChain == 1
		% put h1j into cell arrays
		for i = find(~cellfun('isempty',h1j'))
			H1 = cell(para.nChains,1);
			H1{i} = h1j{i};
			h1jnew = contractMultiChainOBB(BUb, H1, para);
			Hleft = Hleft + updateCleft([],B,[],h1jnew,A,[]);	% add each sub-chain into Hright
		end
	end

	if multiChain == 0
		for m=1:M
            Hleft = Hleft+updateCleft(Opstorage1{m},B,BUb,h2j{m},A,AUb);
		end
	elseif multiChain == 1
		for m=1:M
			h2jnew = contractMultiChainOBB(BUb, h2j(m,:), para);
			Hleft = Hleft+updateCleft(Opstorage1{m},B,[],h2jnew,A,[]);	% transform the interaction terms O_j O_{j+1} in the eff. right basis r_{j-1}
		end
	end

end