function [mps,Vmat,para,op] = tdvp_1site_evolveTTM(mps,Vmat,para,results,op,TT)
%% function tdvp_1site_evolveTTM(mps,Vmat,para,results,op,outFile)
% evolves the last site of the chain with TTM method.
% needs full tmps as input, needs storMPS = 1; read in MPS history from file!
%
% para.timeslice;		% current slice
if ~isfield(op,'BondProjector')
	% No proper preallocation needed, since cell()
	op.BondProjector = cell(para.timeslice,1);		% stores bond-projectors for each past-MPS
end

% Get mps{t-1}
tmps = TT.LastMPS;
tVmat = TT.LastVmat;
	
% Calculate (temporary half-sweep) BondProjector for t-1 -> t
op.BondProjector{para.timeslice} = bondProject(1);

% estimate starting point of TTM
currentN = max(1,para.timeslice-length(TT.TV)+2);	% at max, as long as TTM
rho = 0;
for i = currentN:para.timeslice
	% apply previous full-sweep BondProjector calculated in tdvp_1site, to bring projectors t-2 -> t-1
	if i < para.timeslice-1
		op.BondProjector{i} = op.BondProjector{para.timeslice-1}*op.BondProjector{i};
	end
	% time-evolve past MPS with TTM
	rho = rho + contractRhoTTM(i);						% i is index in op & timeslice
end
% time-evolve current MPS with TTM!
rho = rho + contractRhoTTM(para.timeslice + 1);		%  i=0 indicates MPS

% now need to decompose rho_(l',n~',l,n~)
dr = size(rho);
rho = reshape(rho,dr(1)*dr(2),[]);		%rho_(l',n~'),(l,n~)
[V,D] = eig(rho);
% figure(300); hold all; plot(real(diag(D)),'Displayname',num2str(para.timeslice)); drawnow;	% only debug

% Truncate D and take root
D = real(diag(D));						% rho self-adjoint (hermitian), even positive definite
keepdims = abs(D)>para.svmintol;		% safe threshold?
if sum(keepdims) == 1
	keepdims(2) = 1;
end
V = V(:,keepdims);
D = diag(D(keepdims));			% now: norm(T-V*D*V') < 1e-14
Amat = reshape(V*sqrt(D),dr(1),dr(2),[]);	% A_(l,n,r)
para.D(para.L) = sum(keepdims);	% new right BondDim

% decompose V into mps and Vmat
Amat = permute(Amat,[1,3,2]);	% A_(l,r,n)
dr = size(Amat);
[U,S,V] = svd2(reshape(Amat, [],dr(3)));	% A_(l,r),n = U_(l,r),i * S_i * V_(i,n)
[U, S, V, dr(3), err] = truncateUSV(U, S, V, para, para.d_opt_min);
para.d_opt(para.L) = dr(3);

% only return on-site-MPS
Vmat = V.';						% only transpose! V_(n,n~) = V_(i,n).'
mps = reshape(U*S,dr);

	function rho = contractRhoTTM(t)
		%% function rho = contractRhoTTM(t)
		% contracts the previous mps{t} with the Transfer Tensor to form rho
		% rho(N*dt) is the previous step -> TT{2};
		% rho(t*dt) -> TT{N+2-t}
		% i: index in current Block; t: absolute index
		N = para.timeslice;
		TV = TT.TV{N+2-t};			% TV_((n~',n'),i);  Use copy to save temporary results.
		dL = sqrt(size(TV,1));		% local dimension of site L
		TV = reshape(TV,dL,dL,[]);	% TV_(n~',n',i)
		if t == N+1
			% have time-evolution of current MPS -> no need for Cl or tmps/tVmat as = eye()
			A = mps{end};
			Vm = Vmat{end};
			A = contracttensors(A,3,3,Vm,2,2);
			rho = contracttensors(A,3,2,conj(A),3,2);
			return
		else
			A  = TT.EndMPS{t,1};
			Vm = TT.EndMPS{t,2};
			if t < N
				Cl = op.BondProjector{N}*op.BondProjector{t};		% project from t-1 -> t (half sweep)
			else
				Cl = op.BondProjector{t};							% this is already the current Projector
			end
		end
		
		% 1. contract Vmat into TV
		TV  = contracttensors(TV,3,2,Vm,2,1);		% TV_(n~',i,nOBB) = TV_(n~',n',i) * V_(n',nOBB)
		if i ~= 0
			% 2. contract Bond Projector into MPS
			A   = contracttensors(Cl,2,2,A,3,1);	% A_(l',r,nOBB) = Cl_(l',l) * A_(l,r,nOBB)
		end
		% 3. contract TV and A
		TV  = contracttensors(A,3,3,TV,3,3);				% TV_(l',r,n~',i) = A_(l',r,nOBB) * TV_(n~',i,nOBB)
		TV  = permute(TV,[1,3,2,4]);						% TV_(l',n~',r,i)
		% 4. contract D into TV
		d = size(TV);
% 		rho = contracttensors(TV,4,4,diag(TT.TD{N+1-t}),2,1);		% rho_(l',n~',r,i) = TV_(l',n~',r,i) * D_(i,i)
		rho = reshape(reshape(TV,[],d(4)) * diag(TT.TD{N+2-t}),d);	% This is 50% faster!
		% 5. contract conj(TV) to form V*D*V'
% 		res = contracttensors(rho,4,[3 4], conj(TV),4,[3 4]);		% rho_(l',n~',l,n~) = rho_(l',n~',r,i) * TV*_(l,n~,r,i)
		rho = reshape(reshape(rho,d(1)*d(2),[])*reshape(TV,d(1)*d(2),[])',d([1,2,1,2]));	% rho_(l',n~',l,n~) = rho_(l',r,n~',i) * TV*_(l,r,n~,i)
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

