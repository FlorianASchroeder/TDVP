function X = contractMultiChainOBB(Vmat, mcOp, para)
	%% function contractMultiChainOBB(mcOp)
	% Contracts a multi-chain OBB with a single chain operator stored
	% in mcOp
	% mcOp is empty cell array of length para.nChains
	%  contains 1 operator for subchain n at mcOp{n}
	%  empty spaces indicate eye()

	% Do not accept all-empty mcOp
	if all(cellfun('isempty',mcOp))
		fprintf('ContractMultiChainOBB: Empty mcOp in loop %d at site %d', para.loop, para.sitej);
		X = 0;
		return
	end

% 	mcOp = squeeze(mcOp);		% get rid of singletons
	Vmat = reshape(Vmat,[para.dk(:,para.sitej)',size(Vmat,2)]);
	% contract all empty parts
	ind	 = find(cellfun('isempty',mcOp));		% finds all empty indices
	ind  = reshape(ind, [1,numel(ind)]);		% need row vector!
	X	 = contracttensors(conj(Vmat),para.nChains+1,ind,Vmat,para.nChains+1,ind);		% dk' x d_opt' x dk x d_opt
	ind	 = ~cellfun('isempty',mcOp);				% finds position to contract
	X	 = contracttensors(X,4,[1,3],mcOp{ind},2,[1,2]);
end