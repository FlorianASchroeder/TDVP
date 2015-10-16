function [mps, Vmat, para, results] = tdvp_1site_evolveKn(mps,Vmat,para,results,op,sitej,Cn,Hn)
%% Evolves the non-site center C(n) following Haegeman 2014
%   - Only contains C = exp(i K(n) dt/2) C;
%   - Varies for l->r and l<-r; uses para.sweepto
%   - Hn input as (l*r*n,l*r*n) 2D-array
%
% Created by Florian Schroeder @ Cambridge 20/10/2014

t = para.tdvp.deltaT./2;

if para.tdvp.imagT
	t = -1i*t;
end

switch para.sweepto
    case'r'
        %% Get Dimensions
        % A(n) = A_(l,r^,n), C(n) = C_(cl,cr);  dim(r^) == dim(cl); dim(r)== dim(cr)
        % Where: - A(n) - C(n) - A(n+1) -
        %           |             |
        %           n            n+1
        [BondDimCLeft, BondDimCRight] = size(Cn);                       % BondDimCRight == true right bond dimension of non-site center

		if para.tdvp.expvCustomNow == 0
			%% prepare K(n) from H(n)
			% through contraction with A(t+dt)
			[BondDimALeft, BondDimARight, OBBDim] = size(mps{sitej});       % if no Vmat, then this = dk
			assert(BondDimCLeft == BondDimARight);

			% K(n)_(rl',r',rl,r) = [A*_(l',rl',n') * H(n)_(l',r',n',l,r,n)]_(rl',r',l,r,n) * A_(l,rl,n)
			Hn = reshape(Hn, [BondDimALeft, BondDimCRight, OBBDim, BondDimALeft, BondDimCRight, OBBDim]);
				% Hn has dimensions like focused A(n)
			Kn = contracttensors(conj(mps{sitej}),3,[1 3], Hn,6,[1 3]);                 % K(n)_(rl',r',l,r,n)
			Kn = contracttensors(Kn,5,[3 5], mps{sitej},3,[1 3]);                       % K(n)_(rl',r',r,rl)
			Kn = permute(Kn,[1,2,4,3]);                                                 % K(n)_(rl',r',rl,r)
			Kn = reshape(Kn,[BondDimCLeft*BondDimCRight, BondDimCLeft*BondDimCRight]);  % K(n)_(rl'r',rlr)
		end
        %% Take and apply Matrix exponential
        % C(n,t) = exp(+ i K(n) dt)_(rl'*r',rl*r) * C(n,t+dt)_(rl*r)
		if para.tdvp.expvCustomNow == 0 && size(Kn,1) <= para.tdvp.maxExpMDim
            Cn = expm( 1i .* Kn .* t) * reshape(Cn,[numel(Cn),1]);
			err = 0;
		else
			if para.tdvp.expvCustomNow
				if ~isfield(op,'chain')
					[op] = H_Eff(mps{sitej}, []  , 'CA', op, para);
				else
					op.HleftAV = op.chain(para.currentChain).Hleft;		% star-MPS
					op.h2jAV   = op.chain(para.currentChain).Opleft;
					% make sure H/Opright is correct
				end
				[Cn,err] = expvCustom(1i*t, 'Kn',...
					reshape(Cn,[numel(Cn),1]), para,op);
			else
				if para.tdvp.expvCustomTestAccuracy				% debug
					[op] = H_Eff(mps{sitej}, []  , 'CA', op, para);
					Cn1 = expvCustom(1i*t, 'Kn',...
						reshape(Cn,[numel(Cn),1]), para,op);
				end
				[Cn,err] = expv(1i*t, Kn,...
					reshape(Cn,[numel(Cn),1]),...
					para.tdvp.expvTol, para.tdvp.expvM);
				if para.tdvp.expvCustomTestAccuracyRMS
					disp(rms(Cn-Cn1));		% debug
				end
			end
		end
% 		results.tdvp.expError(para.timeslice,para.expErrorI) = err; para.expErrorI = para.expErrorI+1;
		results.tdvp.expError(para.timeslice,1) = max(results.tdvp.expError(para.timeslice,1),err);

        Cn = reshape(Cn, [BondDimCLeft, BondDimCRight]);
        clear('Kn', 'Hn');

        %% Multiply Ac(n+1) = C(n) * Ar(n+1)
        % set focus on next site A
        % A_(rl,r,n) = C(n)_(rl,l) * A_(l,r,n)
        mps{sitej+1} = contracttensors(Cn,2,2, mps{sitej+1},3,1);
        clear('Cn');
    case 'l'
        %% Get Dimensions
        % A(n+1) = A_(l^,r,n+1), C(n) = C_(cl,cr);  dim(l^) == dim(cr)
        % Where: - A(n) - C(n) - A(n+1) -
        %           |             |
        %           n            n+1
        [BondDimCLeft, BondDimCRight] = size(Cn);
        [BondDimALeft, BondDimARight, OBBDim] = size(mps{sitej+1});         % needed to reshape H(n+1)
        assert(BondDimCRight == BondDimALeft);

		if para.tdvp.expvCustomNow == 0
			%% prepare K(n) from H(n+1)
			% through contraction with A(t+dt)
			% H(n) still persistent from previous sweep
			% K(n)_(l',lr',l,lr) = [A*_(lr',r',n') * H(n+1)_(l',r',n',l,r,n)]_(lr',l',l,r,n) * A_(lr,r,n)
			Hn = reshape(Hn, [BondDimCLeft, BondDimARight, OBBDim, BondDimCLeft, BondDimARight, OBBDim]);
			Kn = contracttensors(conj(mps{sitej+1}),3,[2 3], Hn,6,[2 3]);                   % K(n)_(lr',l',l,r,n)
			Kn = contracttensors(Kn,5,[4 5], mps{sitej+1},3,[2 3]);                         % K(n)_(lr',l',l,lr)
			Kn = permute(Kn,[2,1,3,4]);                                                     % K(n)_(l',lr',l,lr)
			Kn = reshape(Kn,[BondDimCLeft*BondDimCRight, BondDimCLeft*BondDimCRight]);      % K(n)_(l'*lr',l*lr)
		end

        %% Take and apply Matrix exponential
        % C(n,t+dt/2) = exp(+ i K(n) dt/2)_(l'*lr',l*lr) * C(n,t+dt)_(l*lr)
		if para.tdvp.expvCustomNow == 0
			if size(Kn,1) <= para.tdvp.maxExpMDim
				Cn = expm( 1i .* Kn .* para.tdvp.deltaT./2) *reshape(Cn,[numel(Cn),1]);
				err = 0;
			else
				if para.tdvp.expvCustomTestAccuracy								% debug
					para.sitej = para.sitej+1;								% needed for multi-chain reshape
					[op] = H_Eff(mps{sitej+1}, []  , 'CA', op, para);
					Cn1 = expvCustom(1i*t, 'Kn',...
						reshape(Cn,[numel(Cn),1]), para,op);
					para.sitej = para.sitej-1;
				end
				[Cn,err] = expv(1i*t, Kn,...
					reshape(Cn,[numel(Cn),1]),...
					para.tdvp.expvTol, para.tdvp.expvM);
				if para.tdvp.expvCustomTestAccuracyRMS
					disp(rms(Cn-Cn1));		% debug
				end
			end
		else
			para.sitej = para.sitej+1;								% needed for multi-chain reshape
			[op] = H_Eff(mps{sitej+1}, []  , 'CA', op, para);
			[Cn,err] = expvCustom(1i*t, 'Kn',...
				reshape(Cn,[numel(Cn),1]), para,op);
			para.sitej = para.sitej-1;
		end
% 		results.tdvp.expError(para.timeslice,para.expErrorI) = err; para.expErrorI = para.expErrorI+1;
		results.tdvp.expError(para.timeslice,1) = max(results.tdvp.expError(para.timeslice,1),err);

        Cn = reshape(Cn, [BondDimCLeft, BondDimCRight]);
        clear('Kn', 'Hn');

		if sitej == 0
			mps = Cn;		% shortcut for StarMPS
			return;
		end

        %% Multiply Ac(n) = Al(n) * C(n)
        % set focus on next site A
        % A_(l,r,n) = A_(l,r',n) * C(n)_(r',r)
        mps{sitej} = contracttensors(mps{sitej},3,2, Cn,2,1);               % A_(l,n,r)
        mps{sitej} = permute(mps{sitej},[1,3,2]);                           % A_(l,r,n)
        clear('Cn');
end
end