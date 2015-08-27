function [T, dOut] = tensShape(T,select,i,d)
% Use as:
%	tensShape(T, 'unfold', i_unfold, [d1,d2,d2])
%	tensShape(T, 'fold', i_unfold, [d1,d2,d2])
switch lower(select)
	case 'unfold'
		T = tensUnfold(T,i,d);

	case 'fold'
		T = tensFold(T,i,d);
end

	function T = tensUnfold(T,i,d)
	%% T = tensUnfold(T,i)
	%	unfolds tensor by circular permutation
	% dOut gives size(T) after permutation.
		if i ~= 1
			r = length(d);
			T = permute(T,[i:r,1:i-1]);					% permutes only i-1 times
		end
		dOut = size(T);
		T = reshape(T,[d(i),prod(d)/d(i)]);
	end

	function T = tensFold(T,i,d)
	%% T = tensFold(T,i,dim)
	%	folds a tensor which was unfolded along dimension i
	%	d: vector indicates the original dimensions [n1,n2,n3,..];
		r = length(d);
		T = reshape(T,[d(i:r),d(1:i-1)]);
		if i ~= 1
			i = r-(i-1)+1;								% need +1 since [i:r,1:i-1] only shifts by i-1
			T = permute(T,[i:r,1:i-1]);
		end
	end
end