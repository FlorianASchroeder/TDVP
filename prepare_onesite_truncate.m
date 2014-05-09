function [B, U, para, results] = prepare_onesite_truncate(A, direction,para,sitej,results)
% Truncate Bond dimension via SV of A
%
%

[D1, D2, d] = size(A);
D_old=para.D;
switch direction
    case 'lr'
        error('Haven not implemented yet!');
    case 'rl'
        if para.parity=='n'
            A = permute(A, [1, 3, 2]);
            A = reshape(A, [D1, d * D2]);
            [U, S, B] = svd2(A);
            vNE = vonNeumannEntropy(S);
            sv = diag(S);
            if sitej<=min(para.trustsite(end),para.L)
            %Truncate A dims
				keepdims=find(sv>para.svmintol);
				if length(keepdims)>1 && sitej<para.trustsite(end) %%D should be at least 2
						S=S(keepdims,keepdims);
						U=U(:,keepdims);
						B=B(keepdims,:);
						para.D(sitej-1)=length(keepdims);
				end
            end

            DB = size(S, 1);
            B = reshape(B, [DB, d, D2]); B = permute(B, [1, 3, 2]);
            U = U * S;

            if sitej<=para.trustsite(end)
            %Expand A dims
					if sv(end)>para.svmaxtol
						dimincratio=log10(sv(end)/para.svmaxtol)/2;
						adddim=ceil(D1*dimincratio); %increase 20%

						para.D(sitej-1)=para.D(sitej-1)+adddim;
						addmat=zeros(D1,adddim);
						U=cat(2,U,addmat);
						[a,b,c]=size(B);
						addmat=zeros(adddim,b,c);
						B=cat(1,B,addmat);

					end
            end
        else
            error('Haven not implemented different parity yet!');
        end
end

% results.Amat_vNE(sitej) =vNE;
% if length(sv)==para.D || (length(sv)==para.D/2 && (para.parity~='n'))
%     results.Amat_sv(sitej,:)=sv;
% end


end
