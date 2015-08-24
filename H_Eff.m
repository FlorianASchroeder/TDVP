function [op] = H_Eff(mps, V, target, op, para)
%% [op] = H_Eff(mps, Vtarget, op, para)
%	creates the effective operators for 'target' to allow fast contractions
%	in matrix exponentials and eigs()
%
%	mps		: single-site A matrix
%	Vmat	: single-site V matrix
% e.g.
%	[op] = H_Eff(Amat, []  , 'V' , op, para);
%	[op] = H_Eff([]  , Vmat, 'A' , op, para);	% transforms op.h1j, h2j to OBB
%	[op] = H_Eff(Amat, []  , 'CA', op, para);	% relies on previous OBB trafo
%
%   Only MC-TDVP code (HOSVD): 'MC-OBB', 'MC-CV', 'MC-V', 'MC-VS', 'MC-A'
%
%   [op] = H_Eff([]  ,Vtens, 'MC-OBB', op, para);  % transforms all chain h1j/h2j into MCOBB, calls MC-CV
%   [op] = H_Eff([]  ,Vtens, 'MC-CV' , op, para);  % transforms single chain into OBB; defined by para.currentChain
%   [op] = H_Eff([]  ,Vtens, 'MC-V'  , op, para);  % create HnonInt and Hleft/rightAS for single-chain V evolution
%   [op] = H_Eff([]  ,[]   , 'MC-VS' , op, para);  % create h12jAV terms for VS evolution
%   [op] = H_Eff([]  ,Vtens, 'MC-A'  , op, para);  % create h1/2jOBB terms, depends on MC-OBB
%s
%	Created 04/08/2015 by FS
M = para.M;
switch target
	case 'V'
		%% multiply MPS into op.Hright, Hleft, Opright, Opleft
		% HleftA_(n~',n~) = A*_(l',r,n~') [Hl_(l',l) * A_(l,r,n~)]_(l',r,n~)
		op.HleftA = contracttensors(op.Hleft,2,2,mps,3,1);
		op.HleftA = contracttensors(conj(mps),3,[1 2], op.HleftA,3,[1 2]);
		% Hright_(n~',n~) = A*_(l,r',n~') [Hr_(r',r) * A_(l,r,n~)]_(r',l,n~)
		op.HrightA = contracttensors(op.Hright,2,2,mps,3,2);
		op.HrightA = contracttensors(conj(mps),3,[1 2], op.HrightA,3,[2 1]);

		op.OpleftA  = cell(M,1);
		op.OprightA = cell(M,1);

		for k=1:M
			% same contractions as above for Hleft/Hright
			op.OpleftA{k}= contracttensors(op.Opleft{k}, 2,2, mps,3,1);
			op.OpleftA{k}= contracttensors(conj(mps),3,[1 2], op.OpleftA{k},3,[1 2]);

			op.OprightA{k} = contracttensors(op.Opright{k},2,2, mps,3,2);
			op.OprightA{k} = contracttensors(conj(mps),3,[1 2], op.OprightA{k},3,[2 1]);
		end
	case 'A'
		%% multiply Vmat into op.h1j, h2j
		% transform all bare H terms of current sitej into OBB
		% works with multi-chain at OBB level
		if isempty(op.h1j)
			op=gen_sitej_h1h2(op,para,para.sitej);
		end
		if para.nChains == 1
			op.h1jOBB = V' * (op.h1j * V);									% faster and more accurate
			op.h2jOBB = cell(M,2);
			for i=1:M
				op.h2jOBB{i,1} = V' * (op.h2j{i,1} * V);									% faster and more accurate
				op.h2jOBB{i,2} = V' * (op.h2j{i,2} * V);									% faster and more accurate
			end
		else  % nChains > 1
			h1jnew = 0;
			for i = find(~cellfun('isempty',op.h1j'))
				H1 = cell(para.nChains,1);
				H1(i) = op.h1j(i);
				h1jnew = h1jnew + contractMultiChainOBB(V, H1, para);
			end
			op.h1jOBB = h1jnew;

			op.h2jOBB = cell(M,2);
			for i = 1:M
				op.h2jOBB{i,1} = contractMultiChainOBB(V, op.h2j(i,1,:), para);
				op.h2jOBB{i,2} = contractMultiChainOBB(V, op.h2j(i,2,:), para);
			end
		end
	case 'CA'
		%% contract op.Hleft,Hright,Opleft,Opright with op.h1jOBB,h2jOBB and MPS matrix
		% save into op.HleftAV,HrightAV,OprightAV,h2jAV
		switch para.sweepto
			case 'r'	% MPS and V contracted into left operators only
				% HleftAV_(r~',r~) = A*_(l',r~',n~) [Hl_(l',l) * A_(l,r~,n~)]_(l',r~,n~)
				op.HleftAV = updateCleft(op.Hleft, mps, [], [], mps, []);

				% h1jAV_(r~',r~) = A*_(l,r~',n~') [A_(l,r~,n~), h1j_(n~',n~)]_(l,r~,n~')
				op.HleftAV = op.HleftAV + updateCleft([],mps,[],op.h1jOBB,mps,[]);

				% op.Opleft will be summed over and added to HleftAV, since
				% it does not interact across CA
				op.h2jAV	= cell(M,1);
				for k = 1:M
					% OpleftAV_(r~',r~) = A*_(l',r~',n~')[[OPl_(l',l) A_(l,r~,n~)]_(l',r~,n~) h2j_(n~',n~)]_(l',r~,n~')
					op.HleftAV  = op.HleftAV + updateCleft(op.Opleft{k}, mps, [], op.h2jOBB{k,2}, mps, []);

					% h2jAV_(r~',r~) = A*_(l,r~',n~') [A_(l,r~,n~) h2j_(n~',n~)]_(l,r~,n~')
					op.h2jAV{k} = updateCleft([], mps, [], op.h2jOBB{k,1}, mps, []);
				end
			case 'l'	% MPS and V contracted into right operators only
				% HHrightAV_(l~',l~) = A*_(l~',r',n~) [Hr_(r',r) * A_(l~,r,n~)]_(r',l~,n~)
				op.HrightAV = updateCright(op.Hright, mps, [], [], mps, []);

				% h1jAV_(l~',l~) = A*_(l~',r,n~') [A_(l~,r,n~), h1j_(n~',n~)]_(l~,r,n~')
				op.HrightAV = op.HrightAV + updateCright([],mps,[],op.h1jOBB,mps,[]);

% 				op.OprightAV = cell(M,1);		% add to op.HrightAV directly
				op.h2jAV	 = cell(M,1);
				for k = 1:M
					% OprightAV_(l~',l~) = A*_(l~',r',n~')[[A_(l~,r,n~) OPr_(r',r)]_(l~,n~,r') h2j_(n~',n~)]_(l~,r',n~')
					op.HrightAV = op.HrightAV + updateCright(op.Opright{k},mps,[],op.h2jOBB{k,1},mps,[]);

					% h2jAV_(l~',l~) = A*_(l~',r,n~') [A_(l~,r,n~) h2j_(n~',n~)]_(l~,r,n~')
					op.h2jAV{k} = updateCright([],mps,[],op.h2jOBB{k,2},mps,[]);
				end
		end
	case 'MC-V'
		%% contract Tensors of Multi-Chain HOSVD-TDVP for evolution of V{k}
		% builds upon 'V' preparation. -> [op] = H_Eff(Amat, []  , 'V' , op, para);
		% V has to be cell with V{end} = VS the tensor-center, unpermuted: n1 x n2 x ... x nNC x n~
		% Saves in [H/op][left/right]AS
		% leave out ~ for clarity
		NC = para.nChains;
		nc = para.currentChain;                                                                        % # chain I am working on
		nTerms = para.M/NC;                                                                            % # interacting terms per chain
		% 1. contract VS into 4th-order tensor
		IndContract = 1:NC; IndContract(nc) = [];
		VSfullContract = contracttensors(conj(V{end}),NC+1, IndContract, V{end}, NC+1, IndContract);   % VS_(nk', n', nk, n)

		op.HleftAS  = contracttensors(op.HleftA, 2,[1,2], VSfullContract,4,[2,4]);                     % H_(nk',nk)
		op.HrightAS = contracttensors(op.HrightA,2,[1,2], VSfullContract,4,[2,4]);

		% this could collect all terms non-interacting with current chain!
		op.HnonInt  = op.HleftAS + op.HrightAS;
		% h1j terms of other chains
		for mc = 1:NC
			if mc == nc, continue; end;
			IndContract = 1:NC+1; IndContract([mc,nc]) = [];
			VS = contracttensors(conj(V{end}),NC+1, IndContract, V{end}, NC+1, IndContract);
			if mc < nc          % -> VS_(nmc',nk',nmc,nk)
				order = [1,3];
			else                % -> VS_(nk',nmc',nk,nmc)
				order = [2,4];
			end
			op.HnonInt = op.HnonInt + contracttensors(op.h1jMCOBB{mc},2,[1,2], VS,4,order);
		end

		for m = 1:para.M
			mc = ceil(m/nTerms);                        % the chain number for current m
			if nc == mc
				op.OpleftAS{m}  = contracttensors(op.OpleftA{m},  2,[1,2], VSfullContract,4,[2,4]);
				op.OprightAS{m} = contracttensors(op.OprightA{m}, 2,[1,2], VSfullContract,4,[2,4]);
			else % sum into op.HnonInt
				IndContract = 1:NC; IndContract([mc,nc]) = [];
				VS = contracttensors(conj(V{end}),NC+1, IndContract, V{end}, NC+1, IndContract);
				if mc < nc          % -> VS_(nmc',nk',nmc,nk)
					order = [1,3];
				else                % -> VS_(nk',nmc',nk,nmc)
					order = [2,4];
				end
				OpTemp = contracttensors(op.OpleftA{m},  2,[1,2], VS, 6, [3,6]);                       % absorb OpleftA/right since dim: n~ > nk
				op.HnonInt = op.HnonInt + contracttensors(op.h2jMCOBB{m,2,mc},2,[1,2], OpTemp,4,order);
				OpTemp = contracttensors(op.OprightA{m}, 2,[1,2], VS, 6, [3,6]);
				op.HnonInt = op.HnonInt + contracttensors(op.h2jMCOBB{m,1,mc},2,[1,2], OpTemp,4,order);
			end
		end
		% finish with terms: op.HnonInt, op.OpleftAS{m}, op.OprightAS{m}
	case 'MC-CV'
		%% contract to get each single chain h1j/h2j into its OBB.
		% equivalent to creating the H_Eff for VS
		if isempty(op.h1j)
			op = gen_sitej_h1h2(op,para,para.sitej);
		end

		nTerms = para.M/para.nChains;
		nc     = para.currentChain;
		M      = 1:para.M;
		M(ceil(M/nTerms) ~=nc) = [];                         % find which m to operate on
		op.h1jMCOBB{nc} = V{nc}' * (op.h1j{nc} * V{nc});
		for m = M
			op.h2jMCOBB{m,1,nc} = V{nc}' * (op.h2j{m,1,nc} * V{nc});
			op.h2jMCOBB{m,2,nc} = V{nc}' * (op.h2j{m,2,nc} * V{nc});
		end
	case 'MC-OBB'
		%% Create all OBB terms for each chain
		for i = 1:para.nChains
			para.currentChain = i;
			[op] = H_Eff([]  , V, 'MC-CV', op, para);
		end
	case 'MC-VS'
		%% Create the kronecker product terms for faster VS exponential
		d = para.d_opt(:,para.sitej);
		nTerms = para.M/para.nChains;
		for k = 1:para.nChains
			M      = 1:para.M;
			M(ceil(M/nTerms) ~=k) = [];                         % find which m to operate on
			op.h12jAV{k} = kron(eye(d(end)), op.h1jMCOBB{k});
			for m = M
				op.h12jAV{k} = op.h12jAV{k} + kron(op.OpleftA{m} , op.h2jMCOBB{m,2,k});
				op.h12jAV{k} = op.h12jAV{k} + kron(op.OprightA{m}, op.h2jMCOBB{m,1,k});
			end
		end
	case 'MC-A'
		%% Create the OBB terms for A, gathering all chains
		% needs all h1jMCOBB to be up-to-date (MC-OBB)
		% use contractMultiChainOBB with objects from here! Copied from 'A'
		h1jnew = 0;
		for i = find(~cellfun('isempty',op.h1jMCOBB))
			H1 = cell(para.nChains,1);
			H1(i) = op.h1jMCOBB(i);
			h1jnew = h1jnew + contractMultiChainOBB(V{end}, H1, para);
		end
		op.h1jOBB = h1jnew;

		op.h2jOBB = cell(M,2);
		for i = 1:M
			op.h2jOBB{i,1} = contractMultiChainOBB(V{end}, op.h2jMCOBB(i,1,:), para);
			op.h2jOBB{i,2} = contractMultiChainOBB(V{end}, op.h2jMCOBB(i,2,:), para);
		end
end
end