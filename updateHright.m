function [Hright] = updateHright(Hright, h1jOBB, Opstorage2, B, BUb, h2jOBB, A, AUb, M, varargin)
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
% 	multiChain = 0;
% 	if nargin > 9
% 		para = varargin{1};
% 		if para.nChains > 1
% 			multiChain = 1;						% multiChain at OBB level
% 		end
% 	end

% 	if isempty(Opstorage2{1}) 				% initialize matrices if empty
% 		for m=1:M
% 			Opstorage2{m} = zeros(size(B,2));	% zeros(D)? (also for first site); zeors(1) for last
% 		end
% 	end

% 	if multiChain == 0 && any(cellfun('isempty',h2jOBB))					% initialize matrices as 0 if empty (only use for single chains)
% 		for m = find(cellfun('isempty',h2jOBB))
% 			h2jOBB{m} = zeros(size(BUb,1));									% zeros(para.dk(j)) since empty h2j in updateCright -> eye -> wrong! need explicit 0!
% 		end																	% multiChain == 1 does not use updateCright -> leave empty
% 	end

	Hright = updateCright(Hright, B, [], [], A, []);						% transform old Hright from effective basis j+1 into j

	if ~isempty(h1jOBB)
		Hright = Hright + updateCright([],B,[],h1jOBB,A,[]);				% transform H1{j} into effective right basis of j-1
	end

	for m=1:M
		if ~isempty(h2jOBB{m}) && ~isempty(Opstorage2{m})
			Hright = Hright + updateCright(Opstorage2{m},B,[],h2jOBB{m},A,[]);	% transform the interaction terms O_j O_{j+1} in the eff. right basis r_{j-1}
		end
	end
end