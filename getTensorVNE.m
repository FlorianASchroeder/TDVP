function res = getTensorVNE(A)
%% function getTensorVNE(A)
%	calculate the von Neumann entropy for any possible partition of tensor A
%	as a proxy for the entanglement and thus most optimal partition into subtensors.
N = ndims(A);
dA = size(A);
res = {};
k = 1;

% First look at single partitions
% the ones with smallest vNE could be singled out alone
for ii = 1:N
	[B, dB] = tensShape(A,'unfoldiso', ii,dA);
	[~,S,~] = svd2(B);
	res{k,1} = ii;
	res{k,2} = vonNeumannEntropy(S);
	res{k,3} = diag(S);
	fprintf('dim %d: %g\n',res{k,1:2});
	k = k+1;
end

% Second look at pairs
for ii = 1:N
	for jj = ii+1:N
		[B, dB] = tensShape(A,'unfoldiso', [ii,jj],dA);
		[~,S,~] = svd2(B);
		res{k,1} = [ii,jj];
		res{k,2} = vonNeumannEntropy(S);
		res{k,3} = diag(S);
		pairSum = sum([res{res{k,1},2}]);
		fprintf('dim [%d,%d]: %5.3g;\t individual sum: %5.3g; \t rel.Compression: %5.1f\n',res{k,1:2}, pairSum, (1-res{k,2}/pairSum)*100);
		k = k+1;
	end
end

% Third look into triads
for ii = 1:N
	for jj = ii+1:N
		for kk = jj+1:N
			[B, dB] = tensShape(A,'unfoldiso', [ii,jj,kk],dA);
			[~,S,~] = svd2(B);
			res{k,1} = [ii,jj,kk];
			res{k,2} = vonNeumannEntropy(S);
			res{k,3} = diag(S);
			triadSum = sum([res{res{k,1},2}]);
			fprintf('dim [%d,%d,%d]: %5.3g;\t individual sum: %5.3g; \t rel.Compression: %5.1f\n',res{k,1:2}, triadSum, (1-res{k,2}/triadSum)*100);
			k = k+1;
		end
	end
end

% Fourth look into quarts
for ii = 1:N
	for jj = ii+1:N
		for kk = jj+1:N
			for ll = kk+1:N
				[B, dB] = tensShape(A,'unfoldiso', [ii,jj,kk,ll],dA);
				[~,S,~] = svd2(B);
				res{k,1} = [ii,jj,kk,ll];
				res{k,2} = vonNeumannEntropy(S);
				res{k,3} = diag(S);
				triadSum = sum([res{res{k,1},2}]);
				fprintf('dim [%d,%d,%d,%d]: %5.3g;\t individual sum: %5.3g; \t rel.Compression: %5.1f\n',res{k,1:2}, triadSum, (1-res{k,2}/triadSum)*100);
				k = k+1;
			end
		end
	end
end
