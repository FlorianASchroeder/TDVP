function [S, U, sv, para] = hosvd(T, para, dim)
%% [S U sv, para] = hosvd(T, para, dim)
% High Order SVD of a Tensor
%	[S U sv] = HOSVD(T)
%	[S U sv, para] = HOSVD(T, para)
%	[S U sv, para] = HOSVD(T, para, dim)
%
%	T    : Tensor
%	para : VMPS parameters for truncation and dimension tracking: enables truncation
%	dim  : select dimension
%
%	S    : decomposed core tensor
% 	U    : matrices for each dimension, [] for excluded dims
%	sv   : singular values
%
%	eg [S, U, sv, para] = hosvd(randn(5,5,5), para, [0 1 1])


d = size(T);
r = length(d);
if nargin < 2
	useTruncate = false;
end
if nargin < 3
	useTruncate = true;
	dim = ones(1,r);
end
U = cell(1,r);
sv = cell(1,r);
for i = 1:r
	if dim(i)
		A = tensShape(T,'unfold',i,d);
		% SVD in dimension (i)
		[U{i}, sv{i}] = svdTruncate(A);
	else
		U{i} = [];
		sv{i} = [];
	end
end
S = T;
% explicit Tensor multiplications do: 2->3->4...->1
for i = 1:r
	S = tensShape(S,'unfold',2,d);		% circular permute by one each time
	d = d([2:r,1]);						% needed for reshape
	if i == r
		S = U{1}'*S;
	else
		S = U{i+1}'*S;
	end
	d(1) = size(S,1);
	S = reshape(S,d);
end

	function [U,sv] = svdTruncate2(A)
		[U,sv] = svd2(A);
		if ~useTruncate, return; end
		if para.d_opt(i,s) == para.d_opt_min, return; end
		discarddims     = find(sv < para.svmintol);
		if isempty(discarddims), return; end			% no expansion here!

		if sv(discarddims(1)-1) > para.svmaxtol         % if next highest not-discarded element too large (causing expansion in next sweep)
			discarddims = discarddims(2:end);           % always keep one element < svmaxtol
		end
		% if Dif > 0: would remove too many dims; if Dif < 0, no
		% problem so set Dif = 0
		Diff = para.d_opt_min + length(discarddims) - para.d_opt(i,s);	% d_opt - (discard + d_min) = -Diff
		if Diff <= 0      % discarddims does not violate d_opt_min
			Diff = 0;
		end                     % else: discarddims would remove too much by amount = diff;
		dispif('remove dims in d_opt',para.logging)
		para.d_opt(s) = para.d_opt(s)-length(discarddims)+Diff;
		Vmat{s}(:,discarddims(Diff+1:end))=[];                          % only cut the end
		mps{s}(:,:,discarddims(Diff+1:end))=[];
		results.Vmat_sv{s}(discarddims(Diff+1:end))=[];
	end

	function [U,sv] = svdTruncate(A)
		% nicer code than the other!
		[U,sv] = svd2(A);
		sv = diag(sv);
		if ~useTruncate, return; end
		if para.d_opt(i,para.sitej) == para.d_opt_min, return; end
		%% Truncate A dims
        keepdims = find(sv > para.svmintol);
        % If smallest SV too large, keep 1 more
        if (sv(keepdims(end)) > para.svmaxtol)
            if keepdims(end) < length(sv)
                keepdims = [keepdims;keepdims(end)+1];
			else
				return;		% no truncation possible
%                 newBondDim = length(sv);        % = old Dim
                % expand Bond dims!
            end
        end
        if length(keepdims) < para.d_opt_min			% keep at least Dmin bonds
            keepdims = 1:para.d_opt_min;
        end
        if length(keepdims) < length(sv)
            U = U(:,keepdims);                          % keep columns
            sv = sv(keepdims);
			para.d_opt(i,para.sitej) = size(U,2);
        end
	end
end