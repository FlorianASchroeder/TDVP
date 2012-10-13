function [Vmat, V, results,d_opt] = prepare_onesiteVmat(Vmat,para,results,sitej)

[dk, d_opt] = size(Vmat);
if dk>=d_opt
    if para.parity=='n'
    [Vmat, S, V] = svd2(Vmat);
    Vmat_vNE = vonNeumannEntropy(S);
    sv=diag(S);
    d_opt = size(S, 1);
    V=S*V;
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
if nargin==4
    results.Vmat_vNE(sitej) = Vmat_vNE;
    results.Vmat_sv{sitej}=sv;
end
end