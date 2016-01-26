function [mps,Vmat,para] = tdvp_1site_evolveTTM(mps,Vmat,para,results,op,TT,outFile)
%% function tdvp_1site_evolveTTM(mps,Vmat,para,results,op,outFile)
% evolves the last site of the chain with TTM method.
% needs full tmps as input, needs storMPS = 1; read in MPS history from file!
%
% para.timeslice;		% current slice
Cl = cell(para.timeslice,1);				% stores bond-projectors for each past-MPS

% only load max. 100 MB at a time, calculate by taking current mps
vars = whos('mps','Vmat');
varSizeMB = sum([vars.bytes])/1024^2;
nPerBlock = floor(100/varSizeMB);			% n timeslices per block
currentN = 1;								% idx of start of not yet copied
while currentN < para.timeslice
	BlockN = min(nPerBlock, para.timeslice-currentN+1);
	tmps = outFile.tmps(currentN+(0:BlockN-1),:);
	tVmat = outFile.tVmat(currentN+(0:BlockN-1),:);
	
	% contract each slice in outFile with mps and Vmat
	for i = 1:BlockN
		Cl{currentN+i-1} = bondProject(i);		% i as index of tmps, tVmat
		rho = contractRhoTTM(i,currentN+i-1);
	end
	currentN = currentN + BlockN;
end


% only return on-site-MPS
mps = mps{end};
Vmat = Vmat{end};	

	function rho = contractRhoTTM(i,t)
		%% function rho = contractRhoTTM(t)
		% contracts the previous mps{t} with the Transfer Tensor to form rho
		% rho(N*dt) is the current step -> TT{1};
		% rho(t*dt) -> TT{N+1-t}
		% i: index in current Block; t: absolute index
		N = para.timeslice;
		TV = TT.TV{N+1-t};			% TV_((n~',n'),i);  Use copy to save temporary results.
		dL = sqrt(size(TV,1));
		TV = reshape(TV,dL,dL,[]);	% TV_(n~',n',i)
		
		% 1. contract Vmat into TV
		TV  = contracttensors(TV,3,2,tVmat{i,end},2,1);		% TV_(n~',i,nOBB) = TV_(n~',n',i) * V_(n',nOBB)
		% 2. contract Bond Projector into MPS
		A   = contracttensors(Cl{t},2,2,tmps{i,end},3,1);	% A_(l',r,nOBB) = Cl_(l',l) * A_(l,r,nOBB)
		% 3. contract TV and A
		TV  = contracttensors(A,3,3,TV,3,3);				% TV_(l',r,n~',i) = A_(l',r,nOBB) * TV_(n~',i,nOBB)
		% 4. contract D into TV
		d = size(TV);
		rho = contracttensors(TV,4,4,TT.TD{N+1-t},2,1);		% rho_(l',r,n~',i) = TV_(l',r,n~',i) * D_(i,i)
		rho = reshape(reshape(TV,[],d(4)) * TT.TD{N+1-t},d);	% which is faster?
		% 5. contract conj(TV) to form V*D*V'
		rho = contracttensors(rho,4,[2 4], conj(TV),4,[2 4]);	% rho_(l',n~',l,n~) = rho_(l',r,n~',i) * TV*_(l,r,n~,i)
	end

	function Cl = bondProject(t)
		%% function Cl = bondProject(timeslice)
		% Projects the bond D(L-1) of mps(t*dt-1) onto D(L-1) of mps(tnow)
		%
		Cl = [];
		% tmps  = outFile.tmps(t,:);			% read from file
		% tVmat = outFile.tVmat(t,:);
		for kk = 1:para.L-1
			Cl = updateCleft(Cl,mps{kk},Vmat{kk},[],tmps{t,kk},tVmat{t,kk});
		end
		% now Cl is Cl_(l',l) for site L; l': mps(tnow), l: mps(t*dt-1)
	end

end

