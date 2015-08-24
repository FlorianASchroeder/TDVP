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

N = length(mps);
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