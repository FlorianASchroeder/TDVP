function [mps,Vmat,para] = prepare(mps,Vmat,para)
% Does Sweep l->r (left normalization) and r->l (right normalization)
% Calculates either using Vmat or without.
% 1. sweep l->r  (with Vmat); 2.sweep, r->l  (without Vmat)
% for each site:
%	1. prepare_onesiteVmat(Vmat{i}) --> V
%	2. Contract mps{i} and V over spin index
%	3. prepare_onesite(mps) to do SVD on mps A-matrices
%	4. Contract remainder of SVD with mps{i+1}
%
% Changed:
%   FS 20/10/2014: - use para.sweeepto for l or r part.
%	FS 24/08/2015: - Multi-Chain Vmat/Vtens support; MPS normalization added;
%	FS 09/10/2015: - Star-MPS support in prepare_Star_MPS(mps,Vmat,para); Only 'l' sweep on each chain!

N = length(mps);
if para.useStarMPS == 1
	[mps,Vmat,para] = prepare_Star_MPS(mps,Vmat,para);
	return;
elseif para.useTreeMPS == 1
	[mps,Vmat,para] = prepare_Tree_MPS(mps,Vmat,para);
end

para.sweepto = 'r';

for i = 1:N-1
    if para.useVmat == 1
		if para.useVtens == 1 && all(i ~= para.spinposition)
			for k = para.nChains:-1:1
				[Vmat{i}{k}, V] = prepare_onesiteVmat(Vmat{i}{k},para);
				Vmat{i}{end} = contracttensors(V,2,2, Vmat{i}{end}, para.nChains+1, para.nChains);
			end
			Vmat{i}{end}     = reshape(Vmat{i}{end}, [], para.d_opt(end,i));
			[Vmat{i}{end},V] = prepare_onesiteVmat(Vmat{i}{end},para);
			Vmat{i}{end}     = reshape(Vmat{i}{end}, para.d_opt(:,i)');
		else
			[Vmat{i},V] = prepare_onesiteVmat(Vmat{i},para);		% Vmat = U * S * V' ; Vmat := U; V:= S * V'
		end
		mps{i} = contracttensors(mps{i},3,3,V,2,2);                 % = Ai_{l,r,n} * V'_{p,n}; This contraction is defined differently to the paper.
    end
    [mps{i}, U] = prepare_onesite(mps{i},para,i);                   % SVD(Ai_(l,r,n)) = Ai_(l,m,n) * U_(m,r)
    mps{i+1} = contracttensors(U,2,2,mps{i+1},3,1);                 % U_(m,l) * A(i+1)_(l,r,n)
    para=gennonzeroindex(mps,Vmat,para,i);                          % only if parity not 'n'
    para=gennonzeroindex(mps,Vmat,para,i+1);                        % only if parity not 'n'
end
i = N;
if para.useVmat == 1
	if para.useVtens == 1
		for k = para.nChains:-1:1
			[Vmat{i}{k}, V] = prepare_onesiteVmat(Vmat{i}{k},para);
			Vmat{i}{end}    = contracttensors(V,2,2, Vmat{i}{end}, para.nChains+1, para.nChains);
		end
		Vmat{i}{end}     = reshape(Vmat{i}{end}, [], para.d_opt(end,i));
		[Vmat{i}{end},V] = prepare_onesiteVmat(Vmat{i}{end},para);
		Vmat{i}{end}     = reshape(Vmat{i}{end}, para.d_opt(:,i)');
	else
		[Vmat{i},V] = prepare_onesiteVmat(Vmat{i},para);		% Vmat = U * S * V' ; Vmat := U; V:= S * V'
	end
	mps{i} = contracttensors(mps{i},3,3,V,2,2);                 % = Ai_{l,r,n} * V'_{p,n}; This contraction is defined differently to the paper.
end

para.sweepto = 'l';
for i = N:-1:2
	[mps{i}, U] = prepare_onesite(mps{i},para,i);
    mps{i-1}    = contracttensors(mps{i-1}, 3, 2, U, 2, 1);
    mps{i-1}    = permute(mps{i-1}, [1, 3, 2]);
	para = gennonzeroindex(mps,Vmat,para,i);
	para = gennonzeroindex(mps,Vmat,para,i-1);
end
% one more SVD to properly normalise MPS:
[mps{1}, CA, para] = prepare_onesite(mps{1},para,1);
fprintf('Norm was: %g', CA);


end

function [mps,Vmat,para] = prepare_Star_MPS(mps,Vmat,para)
%% normalizes the Star-MPS network on each chain

para.sweepto = 'l';
NC = para.nChains;

for mc = 1:para.nChains
	L = para.chain{mc}.L;
	for ii = L:-1:2
		[Vmat{ii}{mc}, V] = prepare_onesiteVmat(Vmat{ii}{mc},para);		% normalize V
		mps{ii}{mc} = contracttensors(mps{ii}{mc},3,3, V.',2,1);		% focus onto MPS

		[mps{ii}{mc}, U] = prepare_onesite(mps{ii}{mc},para,ii);

		if ii ~= 2
			mps{ii-1}{mc}    = contracttensors(mps{ii-1}{mc}, 3, 2, U, 2, 1);
			mps{ii-1}{mc}    = permute(mps{ii-1}{mc}, [1, 3, 2]);
		else
			mps{1} = contracttensors(mps{1},NC+2,mc+1, U,2,1);
			mps{1} = permute(mps{1}, [1:mc, NC+2, mc+1:NC+1]);
		end
		para.D(mc,ii-1) = size(U,2);
	end

end
	d = size(mps{1});
	[mps{1}, U] = prepare_onesite(reshape(mps{1},[1,prod(para.D(:,1)),para.dk(1,1)]),para,1);
	mps{1} = reshape(mps{1},d);
	fprintf('MPS norm: %g',U);
end

function [treeMPS,Vmat,para] = prepare_Tree_MPS(treeMPS,Vmat,para)
%% function [treeMPS,Vmat,para] = prepare_Tree_MPS(treeMPS,Vmat,para)
%
%	recursively prepare the treeMPS

p = para;					% create working copy of para, to temporarily modify the .use* parameters
pIdx = num2cell(treeMPS.treeIdx+1);

if treeMPS.height == 0
	% this is leaf -> call prepare(mps,Vmat,para)
	p.D = treeMPS.D(2:end-1);
	[treeMPS.mps,treeMPS.Vmat,p] = prepare(treeMPS.mps,treeMPS.Vmat,p);
	para.D{pIdx{:}} = structStarMPS.D;
	para.d_opt{pIdx{:}} = structStarMPS.d_opt;
		
else
	
end
end

