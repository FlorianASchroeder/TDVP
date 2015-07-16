function [Hright] = updateHright(Hright, h1j, Opstorage2, B, BUb, h2j, A, AUb, M, varargin)
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
% Modified
%	FS 16/07/2015: - added multi-chain support in OBB via varargin{1} = para
% called as updateHright(op.Hlrstorage{sitej + 1}, op.h1term(:,sitej), op.Opstorage(:,2,sitej+1),mps{sitej},Vmat{sitej}, op.h2term(:,1,sitej,:), mps{sitej},Vmat{sitej}, M, para);
	if nargin > 9
		para = varargin{1};
		multiChain = 1;						% multiChain at OBB level
	else
		multiChain = 0;
	end

	if isempty(Opstorage2{1}) 				% initialize matrices if empty
		for m=1:M
			Opstorage2{m} = zeros(size(B,2));	% zeros(D)? (also for first site); zeors(1) for last
		end
	end

	if multiChain == 0 && any(cellfun('isempty',h2j))				% initialize matrices as 0 if empty (only use for single chains)
		for m = find(cellfun('isempty',h2j))
			h2j{m} = zeros(size(BUb,1));							% zeros(para.dk(j)) since empty h2j in updateCright -> eye -> wrong! need explicit 0!
		end															% multiChain == 1 does not use updateCright -> leave empty
	end

	Hright = updateCright(Hright, B, BUb, [], A, AUb);				% transform old Hright from effective basis j+1 into j

	if multiChain == 0 && ~isempty(h1j)
		Hright = Hright + updateCright([],B,BUb,h1j,A,AUb);			% transform H1{j} into effective right basis of j-1
	elseif multiChain == 1
		% put h1j into cell arrays
		for i = find(~cellfun('isempty',h1j'))
			H1 = cell(para.nChains,1);
			H1{i} = h1j{i};
			h1jnew = contractMultiChainOBB(BUb, H1, para);
			Hright = Hright + updateCright([],B,[],h1jnew,A,[]);	% add each sub-chain into Hright
		end
	end

	if multiChain == 0
		for m=1:M
			Hright = Hright+updateCright(Opstorage2{m},B,BUb,h2j{m},A,AUb);	% transform the interaction terms O_j O_{j+1} in the eff. right basis r_{j-1}
		end
	elseif multiChain == 1
		for m=1:M
			h2jnew = contractMultiChainOBB(BUb, h2j(m,:), para);
			Hright = Hright+updateCright(Opstorage2{m},B,[],h2jnew,A,[]);	% transform the interaction terms O_j O_{j+1} in the eff. right basis r_{j-1}
		end
	end
end