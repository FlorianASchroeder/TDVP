function para=maxshift(para)
% Calculates maximum possible shift from maximum eigenvalue of x-operator.

% Commented by Florian Schroeder 31/01/2014
% Modified:
% 	FS 31/01/2014:	- added support for folded Chain models
%					- changed loop to start at 2nd site instead of 1st (spin)


para.maxshift=para.shift;
for k=1:para.L
	if k == para.spinposition
		para.maxshift(k) = 0;
	else
		if para.foldedChain == 0
			[bp,bm,n] = bosonop(para.dk(k),0,para.parity);
		elseif para.foldedChain == 1            % In case of Supersite Operators, no need for kron(), as max shift is calculated only for one of the chains.
			[bp,bm,n] = bosonop(sqrt(para.dk(j)),0,para.parity);
		end
		x = sqrt(2)/2.*(bp + bm);
		x = full(x);
		para.maxshift(k) = max(eig(x));
	end
end
end