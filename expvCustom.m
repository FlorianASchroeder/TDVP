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
	case 'HAV'
		AFUN = @HAVmultCV;
	case 'Hn'
		AFUN = @HnmultA;
	case 'Kn'
		AFUN = @KnmultCA;
end

M = para.M;

%%
% anorm, rndoff are estimated using norm(H) later
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
% t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;			% properly approximated using norm(H)
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
  if nstep == 1				% was moved down from initial variable definitions
	  anorm = norm(H);
	  t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;
	  t_step = min( t_out-t_now,t_new );		% re-estimate best t_step
	  rndoff= anorm*eps;
  end
%   fprintf('step: %s, anorm: %g, norm(H): %g, j: %d\n', A, anorm, norm(H), j);		% use to confirm: anorm == norm(H) for the corresponding subspace of v
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
	 elseif ireject == mxrej,
	   error('The requested tolerance is too high.');
	 else
		t_step = gamma * t_step * (t_step*tol/err_loc)^xm;
        s = 10^(floor(log10(t_step))-1);
        t_step = ceil(t_step/s) * s;
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
		[~,newOBBDim] = size(Vmat);
		[~,~,OBBDim]  = size(mps);
% 		CV = reshape(CV,[newOBBDim, OBBDim]);

		% copy OBB into h1j, h2j to match HmultVmat behaviour. This is only
		% calles for Bosons, so is fine!
		op.h1j = op.h1jOBB;
		op.h2j = op.h2jOBB;

		% reuse HmultVmat. Save all temp-results in w to save memory!
		w = HmultVmat(CV, op, newOBBDim,OBBDim, M,para.parity);

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
		[~,OBBDim]		 = size(op.h1jOBB);

		% overhead only for spin site -> negligible!
		op.h1j = op.h1jOBB;
		op.h2j = op.h2jOBB;

		w = HmultA(A, op, BondDimLeft, BondDimRight, OBBDim, M,para.parity,[]);

	end

	function w = KnmultCA(CA)
		%%
		% CA_(r^,r) given as vector, r^ indicates now BondDim from SVD on A
		%	is sweep direction dependent! C_(l,l^) and C_(r^,r)
		% A_(l,r^,n) in mps as tensor.
		[BondDimALeft, BondDimARight, ~] = size(mps);
		[~,BondDimRight] = size(op.Hright);
		[~,BondDimLeft]  = size(op.Hleft);
		w = 0;

		switch para.sweepto
			case 'r'
				CA = reshape(CA, [BondDimARight,BondDimRight]);
				w = w + op.HleftAV * CA + CA * op.Hright.';
				for k = 1:M
					w = w + op.h2jAV{k} * CA * op.Opright{k}.';
				end
			case 'l'
				CA = reshape(CA, [BondDimLeft,BondDimALeft]);
				w = w + op.Hleft * CA + CA * op.HrightAV.';
				for k = 1:M
					w = w + op.Opleft{k} * CA * op.h2jAV{k}.';
				end
		end

		% input V is always vectorized -> reshape
		w = reshape(w, [numel(w),1]);
	end
end

