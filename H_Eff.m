function [op] = H_Eff(mps, V, target, op, para)
%% [op] = H_Eff(target, op, para)
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
				op.HleftAV = contracttensors(op.Hleft,2,2,mps,3,1);
				op.HleftAV = contracttensors(conj(mps),3,[1 3],op.HleftAV,3,[1 3]);

				% h1jAV_(r~',r~) = A*_(l,r~',n~') [A_(l,r~,n~), h1j_(n~',n~)]_(l,r~,n~')
				OpTemp   = contracttensors(mps,3,3, op.h1jOBB.',2,1);
				OpTemp   = contracttensors(conj(mps),3,[1 3], OpTemp,3,[1 3]);
				op.HleftAV = op.HleftAV + OpTemp;

				% op.Opleft will be summed over and added to HleftAV, since
				% it does not interact across CA
				op.h2jAV	= cell(M,1);
				for k = 1:M
					% OpleftAV_(r~',r~) = A*_(l',r~',n~')[[OPl_(l',l) A_(l,r~,n~)]_(l',r~,n~) h2j_(n~',n~)]_(l',r~,n~')
					OpTemp = contracttensors(op.Opleft{k},2,2, mps,3,1);
					OpTemp = contracttensors(OpTemp,3,3, op.h2jOBB{k,2}.',2,1);
					OpTemp = contracttensors(conj(mps),3,[1 3], OpTemp,3,[1 3]);
					op.HleftAV = op.HleftAV + OpTemp;

					% h2jAV_(r~',r~) = A*_(l,r~',n~') [A_(l,r~,n~) h2j_(n~',n~)]_(l,r~,n~')
					op.h2jAV{k}    = contracttensors(mps,3,3, op.h2jOBB{k,1}.',2,1);
					op.h2jAV{k}    = contracttensors(conj(mps),3,[1 3], op.h2jAV{k},3,[1 3]);
				end
			case 'l'	% MPS and V contracted into right operators only
				% HHrightAV_(l~',l~) = A*_(l~',r',n~) [Hr_(r',r) * A_(l~,r,n~)]_(r',l~,n~)
				op.HrightAV = contracttensors(op.Hright,2,2,mps,3,2);
				op.HrightAV = contracttensors(conj(mps),3,[2 3],op.HrightAV,3,[1 3]);

				% h1jAV_(l~',l~) = A*_(l~',r,n~') [A_(l~,r,n~), h1j_(n~',n~)]_(l~,r,n~')
				OpTemp   = contracttensors(mps,3,3, op.h1jOBB.',2,1);
				OpTemp   = contracttensors(conj(mps),3,[2 3], OpTemp,3,[2 3]);
				op.HrightAV = op.HrightAV + OpTemp;

% 				op.OprightAV = cell(M,1);		% add to op.HrightAV directly
				op.h2jAV	 = cell(M,1);
				for k = 1:M
					% OprightAV_(l~',l~) = A*_(l~',r',n~')[[A_(l~,r,n~) OPr_(r',r)]_(l~,n~,r') h2j_(n~',n~)]_(l~,r',n~')
					OpTemp = contracttensors(mps,3,2, op.Opright{k}.',2,1);
					OpTemp = contracttensors(OpTemp,3,2, op.h2jOBB{k,1}.',2,1);
					OpTemp = contracttensors(conj(mps),3,[2 3], OpTemp,3,[2 3]);
					op.HrightAV = op.HrightAV + OpTemp;

					% h2jAV_(l~',l~) = A*_(l~',r,n~') [A_(l~,r,n~) h2j_(n~',n~)]_(l~,r,n~')
					op.h2jAV{k}     = contracttensors(mps,3,3, op.h2jOBB{k,2}.',2,1);
					op.h2jAV{k}     = contracttensors(conj(mps),3,[2 3], op.h2jAV{k},3,[2 3]);
				end
		end


end