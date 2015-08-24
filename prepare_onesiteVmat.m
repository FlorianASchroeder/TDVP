function [Vmat, V, results, d_opt, sv] = prepare_onesiteVmat(Vmat,para,results,sitej,MinDim)
% SVD on Vmat{i} per site. only if dk >= d_opt
% nargin: results, sitej only used for writing output to results.***
% TODO: Use QR to speed up if no result needed
%	put result saving within first if-statement, as otherwise assigned NULL??

[dk, d_opt] = size(Vmat);
if dk>=d_opt
	if para.parity=='n'
		[Vmat, S, V] = svd2(Vmat);				% TODO: if nargin ==4 then SVD, store results else QR; only for 'n'
		sv=diag(S);
		if nargin == 5 && MinDim < d_opt % Do Truncation
			[Vmat, sv, V] = truncateUSV(Vmat, sv, V, para, MinDim);
			S = diag(sv);
		end
        Vmat_vNE = vonNeumannEntropy(S);
        d_opt = size(S, 1);
        V = S*V;
    else
        [Vmat_odd, S_odd, V_odd] = svd2(Vmat(1:dk/2,1:d_opt/2));
        [Vmat_even, S_even, V_even] = svd2(Vmat(dk/2+1:end,d_opt/2+1:end));
        Vmat_vNE_odd = vonNeumannEntropy(S_odd);
        Vmat_vNE_even = vonNeumannEntropy(S_even);
        Vmat_vNE=Vmat_vNE_odd+1i*Vmat_vNE_even;
        Vmat=blkdiag(Vmat_odd,Vmat_even);
        sv_odd=diag(S_odd);
        sv_even=diag(S_even);
        sv=sv_odd+1i*sv_even;
        d_opt = size(S_odd, 1)+size(S_even,1);
        V = blkdiag(S_odd * V_odd,S_even*V_even);
    end
end
if nargin == 4 && ~isempty(results)
    results.Vmat_vNE(sitej) = Vmat_vNE;
    results.Vmat_sv{sitej}  = sv;
end
end