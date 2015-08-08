function [B, U, para, results] = prepare_onesite_truncate(A,para,sitej,results)
% Truncate Bond dimension via SV of A
%
% Changed:
%   FS 20/10/2014: - removed explicit direction, using para.sweepto now!
%   FS 28/10/2014: - implemented sweep to right for TDVP
%                  - Added keeping 1 more dim if smallest SV too large,
%                    limited by para.tdvp.maxBondDim in 'r' sweep
%                  - major restructuring using truncation function
%	FS 09/04/2015: - in 'r' sweep: apply to site L to normalize MPS. No
%					 truncation or para.D saving!

[D1, D2, d] = size(A);
D_old = para.D;
switch para.sweepto
    case 'r'
        if para.parity=='n'
            %% normal Parity
            A = permute(A, [3, 1, 2]);
            A = reshape(A, [d * D1, D2]);		% reshapes to (a1 d),a2
            [B, S, U] = svd2(A);             % Anew with orthonormal columns
            sv = diag(S);
            if para.tdvp.truncateExpandBonds && sitej ~= para.L         % only used in TDVP. For VMPS add: && sitej <= para.trustsite(end)
                if sv(end) > para.svmaxtol
                    %% Expand A dims
                    dimincratio = log10(sv(end)/para.svmaxtol)/2;
                    adddim = ceil(D2*dimincratio);          %increase 20%
                    adddim = min(adddim, para.tdvp.maxBondDim - para.D(sitej));
    %                 para.D(sitej) = para.D(sitej)+adddim;

                    A = cat(2,A,zeros(d * D1,adddim));
                    [B, S, U] = svd2(A);                % better since B contains no zeros but orthonormal vectors!
                    U = U(:,1:D2);                      % form old Bond Dim in columns == U=cat(2,U,addmat);
                else
                    %% Truncate A dims
                    [B, S, U] = truncate(B,S,U);
                end
            end
            DB = size(S, 1);				% new a2 dimension of B
			if sitej ~= para.L
	            para.D(sitej)=DB;
			end
            B = reshape(B, [d, D1, DB]);
            B = permute(B, [2, 3, 1]);
            U = S * U;
        else
            error('Have not implemented different parity yet!');
        end
        resultsSite = sitej;

    case 'l'
        if para.parity=='n'
            A = permute(A, [1, 3, 2]);
            A = reshape(A, [D1, d * D2]);
            [U, S, B] = svd2(A);
            sv = diag(S);
            %% Start Truncation / Expansion
            if sitej <= min(para.trustsite(end),para.L)
                if sv(end) > para.svmaxtol
                    %% Expand A dims
                    % Expand A then do SVD since this leaves B with orthonormal
                    % rows. Padding after SVD would only leave zeros.
                    dimincratio = log10(sv(end)/para.svmaxtol)/2;
                    adddim = ceil(D1*dimincratio);          %increase 20%
%                     para.D(sitej-1) = para.D(sitej-1)+adddim;

                    A = cat(1,A,zeros(adddim, d * D2));
                    [U, S, B] = svd2(A);                % better since B contains no zeros but orthonormal vectors!
                    U = U(1:D1,:);                      % form old Bond Dim in columns == U=cat(2,U,addmat);
                else
                    %% Truncate A dims
                    [U, S, B] = truncate(U,S,B);
                    % If Bond Dim truncated then sv(end) < para.svmaxtol
                end
            end

            DB = size(S, 1);
            para.D(sitej-1) = DB;
            B = reshape(B, [DB, d, D2]); B = permute(B, [1, 3, 2]);
            U = U * S;

        else
            error('Have not implemented different parity yet!');
        end
        resultsSite = sitej-1;
end

if nargin == 4
    results.Amat_vNE(resultsSite) = vonNeumannEntropy(S);
    sv = diag(S);                                               % Expand and truncate might have changed S
    if length(sv)==para.D(resultsSite) || (length(sv)==para.D(resultsSite)/2 && (para.parity~='n'))
        results.Amat_sv{resultsSite} = sv;
    end
end

    function [U, S, V] = truncate(U,S,V)
    %% Truncates Bond Dimensions
    % keeps always one SV in range [para.svmaxtol, para.svmintol]
    %   if not possible -> Do nothing, expand later
    % keeps always at least para.Dmin many dimensions
    % by Florian Schroeder 29/10/2014
        sv = diag(S);
        %% Truncate A dims
        keepdims = find(sv > para.svmintol);
        % If smallest SV too large, keep 1 more
        if (sv(keepdims(end)) > para.svmaxtol)
            if keepdims(end)+1 <= length(sv)
                keepdims = [keepdims;keepdims(end)+1];
            else
%                 newBondDim = length(sv);        % = old Dim
                % expand Bond dims! Happens automatically!
            end
        end
        if length(keepdims) < para.Dmin                     % keep at least Dmin bonds
            keepdims = 1:para.Dmin;
        end
        if length(keepdims) < length(sv)				% && sitej < para.trustsite(end)
            U = U(:,keepdims);                          % keep columns
            S = S(keepdims,keepdims);
            V = V(keepdims,:);                          % keep rows
%             newBondDim = length(keepdims);
        end
    end

end
