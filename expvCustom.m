function [w, err, hump] = expvCustom( t, A, v, mps, Vmat, para, op)
%% A is string telling the operator to construct:
%    A = {'Hn',   'HAA',  'HAV',   'Kn'}
%  v contains the tensor to evolve in time:
%    v = {mps{j}, Vmat{j}, CV{j},   C{j}}
%  mps contains mps A-tensor if v ~= Amat
%    v2= {[],     mps{j},  mps{j},  mps{j}}
%  Vmat contans V-tensor if v ~=Vmat
%    v3= {[],      [],     Vmat{j},  [] }

n = numel(v);
tol = para.tdvp.expvTol;
m   = min(n,para.tdvp.expvM);

%% set Function handle
switch A
	case 'HAA'
		AFUN = @HAAmultV;		% also works with Multi-Chain
		%% non-interacting Hamiltonian terms: Hleft + Hmid + Hright
		% contracted with MPS
		%
		% Hleft_(n~',n~) = A*_(l',r,n~') [Hl_(l',l) * A_(l,r,n~)]_(l',r,n~)
		Hleft = contracttensors(op.Hleft,2,2,mps,3,1);						% Hleft_(l',r,n~) = Hl_(l',l) * A_(l,r,n~)
 		op.HleftA = contracttensors(Hleft,3,[1 2],conj(mps),3,[1 2]);			% Hleft_(n~,n~')  = Hleft_(l',r,n~) * A*_(l',r,n~')
		% Hright_(n~',n~) = A*_(l,r',n~') [Hr_(r',r) * A_(l,r,n~)]_(r',l,n~)
		Hright = contracttensors(op.Hright,2,2, mps,3,2);					% Hright_(r',l,n~) = Hr_(r',r) * A_(l,r,n~)
 		op.HrightA = contracttensors(Hright,3,[2 1], conj(mps),3,[1 2]);		% Hright_(n~,n~')  = Hright_(r',l,n~) * A*_(l',r',n~')
		for k = 1:para.M
			% Opleft_(n~',n~) = A*_(l',r,n~') [Opleft_(l',l) * A_(l,r,n~)]_(l',r,n~)
			Opleft = contracttensors(op.Opleft{k},2,2,mps,3,1);				% Opleft_(l',r,n~) = Opleft_(l',l) * A_(l,r,n~)
			op.OpleftA{k} = contracttensors(Opleft,3,[1 2],conj(mps),3,[1 2]);		% Opleft_(n~,n~')  = Opleft_(l',r,n~) * A*_(l',r,n~')

			% Opright_(n~',n~) = A*_(l,r',n~') [Opright_(r',r) * A_(l,r,n~)]_(r',l,n~)
			Opright = contracttensors(op.Opright{k},2,2, mps,3,2);			% Opright_(r',l,n~) = Opright_(r',r) * A_(l,r,n~)
			op.OprightA{k} = contracttensors(Opright,3,[2 1],conj(mps),3,[1 2]);	% Opright_(n~,n~')  = Opright_(r',l,n~) * A*_(l',r',n~')
		end
	case 'HAV'
		AFUN = @HAVmultCV;
		%% non-interacting Hamiltonian terms: Hleft + Hmid + Hright
		% contracted with MPS
		%
		% Hleft_(n~',n~) = A*_(l',r,n~') [Hl_(l',l) * A_(l,r,n~)]_(l',r,n~)
		Hleft = contracttensors(op.Hleft,2,2,mps,3,1);						% Hleft_(l',r,n~) = Hl_(l',l) * A_(l,r,n~)
 		op.HleftA = contracttensors(Hleft,3,[1 2],conj(mps),3,[1 2]);			% Hleft_(n~,n~')  = Hleft_(l',r,n~) * A*_(l',r,n~')
		% Hright_(n~',n~) = A*_(l,r',n~') [Hr_(r',r) * A_(l,r,n~)]_(r',l,n~)
		Hright = contracttensors(op.Hright,2,2, mps,3,2);					% Hright_(r',l,n~) = Hr_(r',r) * A_(l,r,n~)
 		op.HrightA = contracttensors(Hright,3,[2 1], conj(mps),3,[1 2]);		% Hright_(n~,n~')  = Hright_(r',l,n~) * A*_(l',r',n~')
		for k = 1:para.M
			% Opleft_(n~',n~) = A*_(l',r,n~') [Opleft_(l',l) * A_(l,r,n~)]_(l',r,n~)
			Opleft = contracttensors(op.Opleft{k},2,2,mps,3,1);				% Opleft_(l',r,n~) = Opleft_(l',l) * A_(l,r,n~)
			op.OpleftA{k} = contracttensors(Opleft,3,[1 2],conj(mps),3,[1 2]);		% Opleft_(n~,n~')  = Opleft_(l',r,n~) * A*_(l',r,n~')

			% Opright_(n~',n~) = A*_(l,r',n~') [Opright_(r',r) * A_(l,r,n~)]_(r',l,n~)
			Opright = contracttensors(op.Opright{k},2,2, mps,3,2);			% Opright_(r',l,n~) = Opright_(r',r) * A_(l,r,n~)
			op.OprightA{k} = contracttensors(Opright,3,[2 1],conj(mps),3,[1 2]);	% Opright_(n~,n~')  = Opright_(r',l,n~) * A*_(l',r',n~')
		end
	case 'Hn'
		AFUN = @HnmultA;
	case 'Kn'
		AFUN = @KnmultCA;
		if para.useVmat     % contract H-terms to OBB; also ok for spinsites! since Vmat = eye
			% h1term to OBB, h1j can be rescaled
			% h1j_(n~',n~) = V*_(n',n~') [h1j_(n',n) V_(n,n~)]_(n',n~)
			op.h1j = Vmat' * op.h1j * Vmat;

			% h2term to OBB, h2j can be rescaled
			% h2j_(n~',n~) = V*_(n',n~') [h2j_(n',n) V_(n,n~)]_(n',n~)
			for i=1:size(op.h2j,1);
				op.h2j{i,1} = Vmat' * op.h2j{i,1} * Vmat;
				op.h2j{i,2} = Vmat' * op.h2j{i,2} * Vmat;
			end
		end
end

M = para.M;

%% find largest eigenvalue to estimate norm of functions.
	opts.disp = 0; opts.tol = para.eigs_tol; opts.issym = 1;
	sigma='la';							% largest algebraic eigenvalue
	if para.complex==1
		opts.isreal = false;
		sigma='lr';						% largest real part eigenvalue
	end
	D = eigs(@(x) AFUN(x),n,1,sigma,opts);
% 	fprintf('%.10g\n',diag(D));
% 	w = HAAmultV(v);
% 	return;

%%
anorm = D;						% assume: all operators symmetric -> estimate norm by largest eigenvalue!
rndoff= anorm*eps;
mxrej = 10;  btol  = 1.0e-7;
gamma = 0.9; delta = 1.2;
mb    = m; t_out   = abs(t);
nstep = 0; t_new   = 0;
t_now = 0; s_error = 0;

k1 = 2; xm = 1/m;
w = reshape(v,[n,1]);
normv = norm(w); beta = normv;				% reshape since v = tensor
fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
t_new = t_out-t_now;
t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;			% TODO: est anorm
s = 10^(floor(log10(t_new))-1); t_new = ceil(t_new/s)*s;
sgn = sign(t);

hump = normv;
while t_now < t_out
  nstep = nstep + 1;
  t_step = min( t_out-t_now,t_new );
  V = zeros(n,m+1);
  H = zeros(m+2,m+2);

% Start finding Krylov subspace using v_n = A v_n-1
  V(:,1) = (1/beta)*w;
  for j = 1:m
		p = AFUN(V(:,j));			% AFUN is handle to specific function
%      p = A*V(:,j);
     for i = 1:j					% Modified Gram-Schmidt orthogonalisation
        H(i,j) = V(:,i)'*p;
        p = p-H(i,j)*V(:,i);
     end;
     s = norm(p);
     if s < btol,					% if residual orthogonal vector shorter than btol -> invariant subspace -> step out
        k1 = 0;
        mb = j;
        t_step = t_out-t_now;
        break;
     end;
     H(j+1,j) = s;
     V(:,j+1) = (1/s)*p;
  end;
  if k1 ~= 0,						% i.e. k1 = 2; = 0 if Krylov subspace < m
     H(m+2,m+1) = 1;
     avnorm = norm(AFUN(V(:,m+1)));		% was: norm(A*V(:,m+1));
  end;
  ireject = 0;
  while ireject <= mxrej,			% mxrej = 10 -> 10 tries for good t_step
     mx = mb + k1;							% mb: number of basis vectors -> dim of Krylov
     F = expm(sgn*t_step*H(1:mx,1:mx));
     if k1 == 0,
		err_loc = btol;						% if dim(Krylov) < m -> precise result -> break while
        break;
	 else % estimate local error
        phi1 = abs( beta*F(m+1,1) );
        phi2 = abs( beta*F(m+2,1) * avnorm );
        if phi1 > 10*phi2,
           err_loc = phi2;
           xm = 1/m;
        elseif phi1 > phi2,
           err_loc = (phi1*phi2)/(phi1-phi2);
           xm = 1/m;
        else
           err_loc = phi1;
           xm = 1/(m-1);
        end;
     end;
     if err_loc <= delta * t_step*tol,
        break;
	 else
		t_step = gamma * t_step * (t_step*tol/err_loc)^xm;
        s = 10^(floor(log10(t_step))-1);
        t_step = ceil(t_step/s) * s;
        if ireject == mxrej,
           error('The requested tolerance is too high.');
        end;
        ireject = ireject + 1;
     end;
  end;
  % err_loc small enough @t_step -> apply the exponential
  mx = mb + max( 0,k1-1 );
  w = V(:,1:mx)*(beta*F(1:mx,1));					% calculate resulting w = exp(tA)v
  beta = norm( w );
  hump = max(hump,beta);

  t_now = t_now + t_step;
  t_new = gamma * t_step * (t_step*tol/err_loc)^xm;	% estimate t_step for next loop
  s = 10^(floor(log10(t_new))-1);
  t_new = ceil(t_new/s) * s;

  err_loc = max(err_loc,rndoff);
  s_error = s_error + err_loc;
end;
err = s_error;
hump = hump / normv;

%%	Nested functions
%
%
	function w = HAAmultV(V)
		%% contracts H(n) with A and A* and V
		%   HAA_(n',n~',n,n~) = H(n)_(l',r',n',l,r,n)*A*_(l',r',n~')*A_(l,r,n~)
		%   w(n',n~') = HAA_(n',n~',n,n~) * V_(n,n~)
		% # Contractions: 4+4*M, *: 3 + 4*M

		% input V is always vectorized -> reshape
		[~, ~, OBBDim]  = size(mps);
		[~,dk] = size(op.h1j);
		M = size(op.h2j,1);
		if OBBDim*dk ~= numel(V)
			fprintf('%.10g * %.10g ~= %.10g\n',OBBDim,dk,numel(V));
		end
		V = reshape(V,[dk, OBBDim]);

		w  = V * op.HleftA;							% w_(n',n~')      = 1_(n',n) * V_(n,n~) * Hleft_(n~,n~')
		w  = w + V * op.HrightA;					% w_(n',n~')      = 1_(n',n) * V_(n,n~) * Hright_(n~,n~')
		w  = w + op.h1j * V;											% w_(n',n~') += H1j_(n',n) * V_(n,n~) * 1_(n~,n~')

		for k = 1:M
			w = w + (op.h2j{k,2} * V) * op.OpleftA{k};						% w_(n',n~')      += H2j_(n',n) * V_(n,n~) * Opleft_(n~,n~')
			w = w + (op.h2j{k,1} * V)  * op.OprightA{k};						% w_(n',n~')       += H2j_(n',n) * V_(n,n~) * Opright_(n~,n~')
		end

		w = reshape(w, [numel(w),1]);		% numel faster than a*b
	end

	function w = HAAmultVNew(V)
		%% HAAmultV(V)
		% Needs previous calculation of op.HrightA, op.HleftA, op.OprightA, op.OpleftA
		% works with multi-chain model!
		[~, ~, OBBDim]  = size(mps);
		if iscell(op.h1j)
			dk = prod(cell2mat(cellfun(@(x) size(x,1),op.h1j, 'UniformOutput',false)));		% = dk for multi-chain Hamiltonians
		else
			dk = size(op.h1j,1);
		end

		w = HmultVmat(V, op, dk,OBBDim, M,para.parity);
	end

	function w = HAVmultCV(CV)
		%% For backward evolution of Center V
		% HAV_(n^',n~',n^,n~) = Vmat*_(n',n^')* HAA_(n',n~',n,n~) Vmat_(n,n^)
		%   w(n^',n~') = HAV_(n^',n~',n^,n~) * CV_(n^,n~)
		% expects h1j, h2j in OBB
		% Vmat = Vmat_(n,n^); since focus is taken to CV_(n^,n~)
		[dk,newOBBDim] = size(Vmat);
		[~,~,OBBDim]  = size(mps);
		CV = reshape(CV,[newOBBDim, OBBDim]);

		% reuse HAAmultV. Save all temp-results in w to save memory!
		w = Vmat * CV;									% VmatCV_(n,n~) = Vmat_(n,n^) * CV_(n^,n~)
		w = reshape(w,[numel(w),1]);					% vectorize for use in HAAmultV
		w = HAAmultV(w);								% HAAV_(n',n~') = HAA_(n',n~',n,n~) * VmatCV_(n,n~)
		w = reshape(w,[dk,OBBDim]);						% HAAV_(n',n~')
		w = Vmat' * w;									% w_(n^',n~') = Vmat*_(n',n^') * HAAV_(n',n~') = HAV * CV

		% vectorize output
		w = reshape(w,[numel(w),1]);

	end

	function w = HnmultA(A)
		%% Multiplies H(n) with MPS-A
		%	uses old function HmultA.m since equivalent! Only now for
		%	parity = 'n'
		%	expect op.h1j, op.h2j transformed to OBB if para.useVmat = 1
		%   since para should not be passed to here.

		% input A is always vectorized -> reshape
		[~,BondDimRight] = size(op.Hright);
		[~,BondDimLeft]  = size(op.Hleft);
		[~,OBBDim]		 = size(op.h1j);

		w = HmultA(A, op, BondDimLeft, BondDimRight, OBBDim, M,para.parity,[]);

	end

	function w = KnmultCA(CA)
		%%
		% C_(r^,r) given as vector, r^ indicates now BondDim from SVD on A
		%	is sweep direction dependent! C_(l,l^) and C_(r^,r)
		% A_(l,r^,n) in mps as tensor.
		[BondDimALeft, BondDimARight, OBBDim] = size(mps);
		[~,BondDimRight] = size(op.Hright);
		[~,BondDimLeft]  = size(op.Hleft);
		w = 0;

		switch para.sweepto
			case 'r'
				CA = reshape(CA, [BondDimARight,BondDimRight]);
				CA = contracttensors(mps,3,2,CA,2,1);						% = A_(l,r^,n) * C_(r^,r) = ()_(l,n,r)
				CA = permute(CA, [1,3,2]);
				CA = reshape(CA, [numel(CA),1]);
				CA = HmultA(CA, op, BondDimALeft, BondDimRight, OBBDim, M,para.parity,[]);	% = ()_(l' * r' * n')
% 				CA = HnmultA(CA);
				CA = reshape(CA, [BondDimALeft, BondDimRight,OBBDim]);	% = ()_(l',r',n')
				w = contracttensors(conj(mps),3,[1,3],CA,3,[1,3]);		% = A*_(l',r^',n') * ()_(l',r',n') = C_(r^',r')
			case 'l'
				CA = reshape(CA, [BondDimLeft,BondDimALeft]);
				CA = contracttensors(CA,2,2,mps,3,1);						% = C_(l,l^) * A_(l^,r,n) = ()_(l,r,n)
				CA = reshape(CA, [numel(CA),1]);
				assert(BondDimARight == BondDimRight);
				CA = HmultA(CA, op, BondDimLeft, BondDimARight, OBBDim, M,para.parity,[]);	% = ()_(l' * r' * n')
% 				CA = HnmultA(CA);											% = ()_(l' * r' * n')
				CA = reshape(CA, [BondDimLeft, BondDimARight,OBBDim]);	% = ()_(l',r',n')
				w = contracttensors(CA,3,[2,3],conj(mps),3,[2,3]);		% = ()_(l',r',n') * A*_(l^',r',n') = C_(l',l^')
		end

		% input V is always vectorized -> reshape
		w = reshape(w, [numel(w),1]);
	end
end

