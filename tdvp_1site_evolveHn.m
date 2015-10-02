function [mps, Vmat, para, results, op, Hn] = tdvp_1site_evolveHn(mps,Vmat,para,results,op,sitej)
%% Evolves one site following Haegeman 2014
%   - Only contains A = exp(-i H(n) dt/2) A;
%   - Same procedure for l->r or l<-r
%   - Splits time-evolution into A and V part if para.useVmat == 1
%
% Created by Florian Schroeder @ Cambridge 20/10/2014
%
% Changed:
%	- FS 22/11/2014: - replaced extra h1j and h2j by op.h1j, op.h2j

[BondDimLeft, BondDimRight, OBBDim]  = size(mps{sitej});
dk = prod(results.dk{end}(:,sitej));
if ~para.useVmat || any(sitej == para.spinposition)
    assert(OBBDim == dk) ;
	if para.nChains > 1
		op = H_Eff([]  , Vmat{sitej}, 'A' , op, para);	% deals with Multi-Chain magic for spinsites
	else
		op.h1jOBB = op.h1j;								% Overhead only for spinsites mostly -> negligible!
		op.h2jOBB = op.h2j;
	end
end
M = size(op.h2j,1);

if sitej ~= para.L
    t = para.tdvp.deltaT./2;
else
    t = para.tdvp.deltaT;
end

if para.tdvp.imagT
	t = -1i*t;
end

%% If using Vmat, evolve it first, only BOSON!
if para.useVmat == 1 && prod(sitej ~= para.spinposition)                % if bosonic site only!
    %% expand OBB in A and V by 50%
    % since this expansion is temporarily, save change later.
    % expand always BEFORE SVD
    if (dk > OBBDim) && para.tdvp.expandOBB
        % next line: argument ,BondDimLeft*BondDimRight-OBBDim in min() is
        % wrong I think. can be removed, but has to be checked again!
		expandBy = min([floor(OBBDim*0.2),para.tdvp.maxOBBDim-OBBDim,BondDimLeft*BondDimRight-OBBDim, dk-OBBDim]);
		if expandBy == 0, expandBy = 1; end
        mps{sitej} = cat(3,mps{sitej},zeros(BondDimLeft, BondDimRight, expandBy));
        Vmat{sitej} = cat(2,Vmat{sitej}, zeros(dk, expandBy));
        [~, ~, OBBDim]  = size(mps{sitej});
        para.d_opt(sitej) = OBBDim;
		para.d_optnew(sitej) = OBBDim;
%     else
%         Anew = mps{sitej};
%         Vmatnew = Vmat{sitej};
    end

    %% SVD to set focus on Vmat
    [Amat,V] = prepare_onesiteAmat(mps{sitej},para,sitej);              % left-normalize A, SVD in n.
    [BondDimLeft, BondDimRight, OBBDim]  = size(Amat);
%	Vmat_focused = Vmat{sitej} * transpose(V);							% set focus on Vmat: V_(n,n~)
	Vmat_focused = contracttensors(Vmat{sitej}, 2, 2, V, 2, 2);

    % Amat = MPS{sitej} left normalised;

	%% if HAA or any other operators have size < 1GB, construct them explicitly.
	% 16 bytes per complex number -> 2^26 = 6.7e7 elements in matrix.
	if (BondDimLeft * BondDimRight * OBBDim) > para.tdvp.maxExpVDim && para.tdvp.expvCustom
		% Largest Operator for matrix Exp. is Hn -> use as criterion
		% use expvCustom for this entire site sweep!
		para.tdvp.expvCustomNow = 1;
	else
		para.tdvp.expvCustomNow = 0;
	end

	op = H_Eff(Amat, []  , 'V' , op, para);		% create effective H terms for V

	if para.tdvp.expvCustomNow == 0
		% now: Construct
		% HAA_(n',n~',n,n~) = H(n)_(l',r',n',l,r,n)*A*_(l',r',n~')*A_(l,r,n~)
		% for: V(t+dt) = exp(-i HAA dt)_(n',n~',n,n~) * V(t)_(n,n~)

		HAA = 0;    % kron( zeros(OBBDim), zeros(dk))

		%% non-interacting Hamiltonian terms: Hleft + Hmid + Hright
		% contracted with MPS
		%
		% HleftA_(n~',n~) = A*_(l',r,n~') [Hl_(l',l) * A_(l,r,n~)]_(l',r,n~)
		HAA = HAA + kron(op.HleftA, eye(dk));

		% HrightA_(n~',n~) = A*_(l,r',n~') [Hr_(r',r) * A_(l,r,n~)]_(r',l,n~)
		HAA = HAA + kron(op.HrightA, eye(dk));

		% Hmid_(n~',n~) = A*_(l,r,n~') * A_(l,r,n~)
		Hmid = contracttensors(conj(Amat),3,[1 2], Amat,3,[1 2]);           % should be eye(OBBDim) TODO: could be replaced.
		HAA = HAA + kron(Hmid, op.h1j);
		clear('Hmid');

		%% Interacting Hamiltonian terms: \sum_i^M op.Opleft
		for m = 1:M
			% Opleft_(n~',n~) = A*_(l',r,n~') [Opleft_(l',l) * A_(l,r,n~)]_(l',r,n~)
			HAA = HAA + kron(op.OpleftA{m}, op.h2j{m,2});
			% Opright_(n~',n~) = A*_(l,r',n~') [Opright_(r',r) * A_(l,r,n~)]_(r',l,n~)
			HAA = HAA + kron(op.OprightA{m}, op.h2j{m,1});
		end
	end
    %% Take matrix exponential
    % V(t+dt) = exp(-i HAA dt)_(n',n~',n,n~) * V(t)_(n,n~)
	if para.tdvp.expvCustomNow == 0
		if size(HAA,1) <= para.tdvp.maxExpMDim
			Vmat_focused = expm(- 1i .* HAA .*t) * reshape(Vmat_focused,[dk*OBBDim,1]);
			err = 0;
		else
			% Do approximation of exp(A)*v
			if para.tdvp.expvCustomTestAccuracy
				V1 = expvCustom(- 1i*t,'HAA',...
					 reshape(Vmat_focused,[dk*OBBDim,1]), para, op);
			end
			[Vmat_focused,err] = expv(- 1i*t,HAA,...
						   reshape(Vmat_focused,[dk*OBBDim,1]),...
						   para.tdvp.expvTol, para.tdvp.expvM);
 			if para.tdvp.expvCustomTestAccuracyRMS
				disp(rms(Vmat_focused-V1));
			end
		end
	else
		[Vmat_focused, err] = expvCustom(- 1i*t,'HAA',...
					   reshape(Vmat_focused,[dk*OBBDim,1]), para, op);
	end
% 	results.tdvp.expError(para.timeslice,para.expErrorI) = err; para.expErrorI = para.expErrorI+1;
	results.tdvp.expError(para.timeslice,1) = max(results.tdvp.expError(para.timeslice,1),err);
    Vmat_focused = reshape(Vmat_focused,[dk,OBBDim]);
%     clear('HAA');
%% TODO: Introduce decay for last site!


    %% normalise Vmat and take focus to A
    [Vmat{sitej}, V, results] = prepare_onesiteVmat(Vmat_focused,para,results,sitej);
	if para.tdvp.expandOBB
		% remove empty SV:
		keep = results.Vmat_sv{sitej} ~= 0;
		results.Vmat_sv{sitej} = results.Vmat_sv{sitej}(keep);
		Vmat{sitej} = Vmat{sitej}(:,keep); V = V(keep, :);
	end
	[n1, n2] = size(V);
	OBBDimNew = n1;
	% put h1j, h2j into OBB. writes into op.h1jOBB, op.h2jOBB only! Since h1j, h2j should be
	op = H_Eff([]  , Vmat{sitej}, 'A' , op, para);

    % V_(n^,n~)
    % evolve center backward in time:
    % HAV_(n^',n~',n^,n~) = Vmat*_(n',n^')* HAA_(n',n~',n,n~) Vmat_(n,n^)

	if para.tdvp.expvCustomNow == 0
		%%
		HAA = reshape(HAA,[dk,OBBDim,dk,OBBDim]);
		HAV = contracttensors(conj(Vmat{sitej}),2,1, HAA,4,1);
		HAV = contracttensors(HAV,4,3, Vmat{sitej},2,1);
		HAV = permute(HAV,[1,2,4,3]);
		[n1,n2,n3,n4] = size(HAV);
		HAV = reshape(HAV, [n1*n2,n3*n4]);
		if size(HAV,1) <= para.tdvp.maxExpMDim
% 			V = expm( 1i.* para.tdvp.deltaT./2.*HAV) * reshape(V,[numel(V),1]);
			[V,err] = expv(+ 1i*t,HAV,...
					reshape(V,[numel(V),1]),...
					para.tdvp.expvTol, para.tdvp.expvM);
		else
			if para.tdvp.expvCustomTestAccuracy
				V1 = expvCustom(+ 1i*t,'HAV',...
					reshape(V,[numel(V),1]), para,op);
			end
			[V,err] = expv(+ 1i*t,HAV,...
					reshape(V,[numel(V),1]),...
					para.tdvp.expvTol, para.tdvp.expvM);
			if para.tdvp.expvCustomTestAccuracyRMS
				disp(rms(V-V1));
			end
		end
	else
		[V,err] = expvCustom(+ 1i*t,'HAV',...
				reshape(V,[numel(V),1]), para,op);
	end
% 	results.tdvp.expError(para.timeslice,para.expErrorI) = err; para.expErrorI = para.expErrorI+1;
	results.tdvp.expError(para.timeslice,1) = max(results.tdvp.expError(para.timeslice,1),err);
    V = reshape(V,[n1,n2]);
    mps{sitej} = contracttensors(Amat, 3, 3, V, 2, 2);
    clear('Amat','Vmat_focused','V');
	OBBDim = OBBDimNew;
	para.d_opt(sitej) = OBBDim;

end

%% Now: construct H(n)
% according to Haegeman 2014
% For Hamiltonian use h1jOBB, h2jOBB

if para.tdvp.expvCustomNow == 0
	if para.tdvp.expvCustomTestAccuracy
		tempT = tic;
	end

% Construct Hn explicitly
	Hn=0;       % Hn = kron(eye(OBBDim),kron(eye(BondDimRight),eye(BondDimLeft)))
	% all terms:
	Hn = Hn + kron(eye(OBBDim),kron(eye(BondDimRight),op.Hleft));
	Hn = Hn + kron(eye(OBBDim),kron(op.Hright,eye(BondDimLeft)));
	Hn = Hn + kron(op.h1jOBB,kron(eye(BondDimRight),eye(BondDimLeft)));
	for m=1:M
		Hn = Hn + kron(op.h2jOBB{m,2},kron(eye(BondDimRight),op.Opleft{m}));
		Hn = Hn + kron(op.h2jOBB{m,1},kron(op.Opright{m},eye(BondDimLeft)));
	end
	if para.tdvp.expvCustomTestAccuracy
		t4 = toc(tempT);
	end
end
%% Take and apply Matrix exponential
% A(t+dt) = exp(-i Hn dt)_(l',r',n',l,r,n) * A(t)_(l,r,n)
% Last site special, see Haegeman 2014
% TODO: change expm() with threshold

if  para.tdvp.expvCustomNow == 0
	if size(Hn,1) <= para.tdvp.maxExpMDim
		mpsNew = expm(- 1i .* Hn .*t) * reshape(mps{sitej},[numel(mps{sitej}),1]);
		err = 0;
	else
		if para.tdvp.expvCustomTestAccuracy									% debug
			tempT = tic;
			mpsNew1 = expm(- 1i .* Hn .*t) * reshape(mps{sitej},[numel(mps{sitej}),1]);
			t1 = toc(tempT);

			tempT = tic;
			mpsNew2 = expvCustom(- 1i*t, 'Hn',...
					  reshape(mps{sitej},[numel(mps{sitej}),1]), para,op);
			t3 = toc(tempT);
		end
		tempT = tic;
		[mpsNew,err] = expv(- 1i*t, Hn,...
				 reshape(mps{sitej},[numel(mps{sitej}),1]),...
				 para.tdvp.expvTol, para.tdvp.expvM);
		t2 = toc(tempT);
		if para.tdvp.expvCustomTestAccuracyRMS
			disp(rms(mpsNew-mpsNew1));	% debug
		end
		if para.tdvp.expvCustomTestAccuracy
			results.tdvp.expvTime = [results.tdvp.expvTime; t1, t2, t3, t4, numel(Hn)]; % times for [ expM, expV, exvCustom, Hn building]
		end
	end
else
	tempT = tic;
	[mpsNew,err] = expvCustom(- 1i*t, 'Hn',...
		reshape(mps{sitej},[numel(mps{sitej}),1]), para,op);
	Hn = [];		% dummy return value;
	t1 = toc(tempT);
	if para.tdvp.expvCustomTestAccuracy
		results.tdvp.expvTime = [results.tdvp.expvTime; t1,0,0,0, BondDimLeft*BondDimRight*OBBDim];
	end
end
% results.tdvp.expError(para.timeslice,para.expErrorI) = err; para.expErrorI = para.expErrorI+1;
results.tdvp.expError(para.timeslice,1) = max(results.tdvp.expError(para.timeslice,1),err);

mps{sitej} = reshape(mpsNew,[BondDimLeft,BondDimRight,OBBDim]);
% now: A and V are time-evolved, A is focused
% if sitej = L, then start lr sweep with decomposition of mps

% Only return current-site matrices!
mps = mps{sitej};
Vmat = Vmat{sitej};


end