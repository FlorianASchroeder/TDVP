function out = getObservable(type,mps,Vmat,para)
%% Calculates the following Observables:
%	SBM:    'spin', 'occupation', 'current', 'shift', 'rdm', 'staroccupation', 'starpolaron',
%			'energy', 'bath1correlators', 'bath2correlators'
%   MLSBM:  'participation', 'tunnelenergy', 'staroccupation', 'starpolaron', 'energy'
%
%   Use as: first argument is cell, with descriptor as type{1}
%               all other type{n} are options defined in each function
%           getObservable({'spin'}					,mps,Vmat,para)
%           getObservable({'occupation'}			,mps,Vmat,para)
%           getObservable({'current',AnAm}			,mps,Vmat,para)
%           getObservable({'shift'}					,mps,Vmat,para)
%           getObservable({'rdm',2}					,mps,Vmat,para)
%           getObservable({'tunnelenergy',op}		,mps,Vmat,para)
%           getObservable({'bath1correlators'}		,mps,Vmat,para)
%           getObservable({'bath2correlators'}		,mps,Vmat,para)
%           getObservable({'starpolaron'}			,mps,Vmat,para)
%           getObservable({'staroccupation',AnAm}	,mps,Vmat,para)
%           getObservable({'energy',op}				,mps,Vmat,para)
%           getObservable({'hshi',op}				,mps,Vmat,para)		% op is optional
%
% In 'current' and 'staroccupation': AnAm optional!
%
%   The output 'out' can be of different form:
%       'spin':             struct, fields: sx, sy, sz
%       'occupation':       array	(L x NC)
%       'shift':            array
%       'rdm':              matrix
%       'participation':    number
%       'tunnelenergy':     number
%		'bath1correlators':	vector	(L x 1)
%		'bath2correlators':	matrix	(L x L x NC)
%		'staroccupation':	array	(3 x X x NC)
%		'current':			array	(NC x L-1)
%		'starpolaron':		array	(3 x X x NC)
%		'energy'			scalar
%		'hshi'				scalar
%		'coherence'			not implemented yet
%
%   Created 03/06/2014 by Florian Schroeder @ University of Cambridge
%   TODO:   - Implement Boson Site Shift to export the routine from optimizesite.m. Only get a single shift value.
%
% Modified:
%	- 21/12/14 FS: replaced OBB contractions with faster matrix products.
%	- 21/08/15 FS: added Multi-Chain Vmat / Vtens capability
switch type{1}
    case 'spin'
        % applicable for spin-boson model and for folded SBM2
        out = calSpin(mps,Vmat,para);

    case 'occupation'
        % applicable to all single stranded chains. Extensible to 2-chains
		if para.foldedChain == 0
			if para.nChains == 1
				out = calBosonOcc(mps,Vmat,para);
			else										% Multi-chain with Vmat / V-tensor-network
				out = calBosonOcc_MC(mps,Vmat,para);	% (NC x L)
			end
		elseif para.foldedChain == 1
			out(1,:) = calBosonOcc(mps,Vmat,para,1);
			out(2,:) = calBosonOcc(mps,Vmat,para,2);
		end
		if size(out,1) > 1
			out = out.';			% (L x NC) for compatibility with tresults
		end

    case 'shift'
        % applicable to all single stranded chains. Extendable to 2-chains
        out = calBosonShift(mps,Vmat,para);

	case '1siteshift'
		% intended for use in optimizesite. Only for focused single-site
		% mps and Vmat are single-site matrices
		out = cal1siteShift(mps,Vmat,para);
		
    case 'rdm'
        % applicable to all systems.
		if type{2} <= para.L
            out = calRDM(mps,Vmat,para,type{2});
		else
            out = [];
		end
		
	case 'rdm_adiabatic'
		% applicable to site 1 for now only!
		if type{2} ~= 1
			error('rdm_adiabatic only available for site 1');
		end
		if para.nChains > 1
			error('rdm_adiabatic only for single chains!');
		end
		assert( type{3} <= para.dk(1,1), 'Need valid number for state selection');
		k = type{3};		% the (adiabatic) bond state to construct the (adiabatic) RDM from
		
		% SVD of mps{1} and delete unwanted bond states!
		A = mps{1};
		[D1, D2, d] = size(A);				% d is site dimension
		A = permute(A, [3, 1, 2]);
		A = reshape(A, [d * D1, D2]);		% reshapes to (a1 d),a2
		[B, S, U] = svd2(A);				% Could also use QR decomposition if nargin !=5
		DB = size(S,1);
		B = B* sparse(k,k,1,DB,DB) *U;
		B = reshape(B, [d, D1, D2]);
		mps{1} = permute(B, [2, 3, 1]);		% now only contains the desired bond state
		
		out = calRDM(mps,Vmat,para,type{2});
		
	case 'rdm_adiabatic2'
		% applicable to site 1 for now only!
		% A different way of selecting the bond-state, without SVD.
		% Seems to be equivalent to 'rdm_adiabatic'
		if type{2} ~= 1
			error('rdm_adiabatic only available for site 1');
		end
		if para.nChains > 1
			error('rdm_adiabatic only for single chains!');
		end
		assert( type{3} <= para.dk(1,1), 'Need valid number for state selection');
		k = type{3};		% the (adiabatic) bond state to construct the (adiabatic) RDM from
		
		[~, D2, ~] = size(mps{1});				% d is site dimension
		for ii = 1:D2
			if ii ~= k
				mps{1}(1,ii,:) = 0;				% delete 2 for k=1 or 1 for k=2
			end
		end
		mps{1} = mps{1}./norm(reshape(mps{1},[],1));		% normalise the state!
		
		out = calRDM(mps,Vmat,para,type{2});
		
    case 'participation'
        % applicable to all systems, only for first site.
        out = calParticipation(calRDM(mps,Vmat,para,1));
		
    case 'tunnelenergy'
        % useful for MLSBM, applicable to all
        % needs type{2} = op
		if length(type) == 2
            out = calTunnelingEnergy(mps,Vmat,para,type{2});
		else
            out = calTunnelingEnergy(mps,Vmat,para);
		end

	case 'bath2correlators'
		% needed for mapping from chain to star
		% returns triangular L x L Matrix
		if para.useStarMPS
			NC = para.nChains;
			out = zeros(para.L, para.L, NC);
			for mc = 1:NC
				% extract MPS for each chain and calculate <anam> separately
				cL = para.chain{mc}.L;
				mpsC  = [mps{1}, cellfun(@(x) x{mc},mps(2:cL),'UniformOutput',false)];
				VmatC = [Vmat{1}, cellfun(@(x) x{mc},Vmat(2:cL),'UniformOutput',false)];
				mpsC{1} = permute(mpsC{1}, [1:mc,mc+2:NC+1, mc+1, NC+2]);
				mpsC{1} = reshape(mpsC{1}, [],para.D(mc,1),para.dk(1,1));

				paraC = para;
				paraC.nChains = 1; paraC.L = para.chain{mc}.L;
				paraC.dk = para.dk(mc,:); paraC.dk(1) = para.dk(1,1);	% need to preserve system dk
				paraC.shift = para.shift(mc,:);
				% calculate <anam>
				out(1:paraC.L,1:paraC.L,mc) = calBath2SiteCorrelators_MC(mpsC,VmatC,paraC);
			end
		else
			out = calBath2SiteCorrelators_MC(mps,Vmat,para);
		end

	case 'staroccupation'
		%% does mapping chain -> star
		% uses bath2correlators
		if para.foldedChain == 1, error('Not yet implemented'); end
		% Get the correlator
		if length(type) == 2
			AmAn = real(type{2});
		else
			AmAn = real(getObservable({'bath2correlators'},mps,Vmat,para));		% L x L x nChains
		end
		% preallocate output array
		[x,~,~] = getStarMapping(para,1);
		out      = zeros(2,length(x),para.nChains);
		for mc = 1:para.nChains
			[x,hsquared,pxn] = getStarMapping(para,mc);
			if strcmp(para.chain{mc}.mapping, 'LanczosTriDiag')
				result = diag(pxn * (AmAn(:,:,mc) + triu(AmAn(:,:,mc),1).') * pxn.');	% pxn == U map -> no need for hsquared
				out(:,1:length(x),mc) = [x,result(2:end)].';
			else
% 				out(:,1:length(x),mc) = [x, hsquared.*((pxn.^2)*diag(AmAn(:,:,mc)) + 2.* diag(pxn*(AmAn(:,:,mc)-diag(diag(AmAn(:,:,mc))))*pxn.') )]';
				out(:,1:length(x),mc) = [x, hsquared.*diag(pxn*(AmAn(:,:,mc)+triu(AmAn(:,:,mc),1).')*pxn.')].';
			end
		end

	case 'current'
		%% gets current through each bond.
		% works with Single and Multi-Chain Vmat / Vtens
		L = para.L; NC = para.nChains;
		if length(type) == 2
			AmAn = imag(type{2});												% (L x L x NC)
		else
			AmAn = imag(getObservable({'bath2correlators'},mps,Vmat,para));		% tridiagonal
		end
		out = zeros(L-1,NC);
		for mc = 1:NC
			out(:,mc) = (para.chain{mc}.t.*diag(AmAn(1:end-1,2:end,mc)));		% t(n)*a(n)*a(n+1)^+
		end

	case 'bath1correlators'
		% = state projected chain polaron
		% needed for mapping from chain to star, starpolaron
		% returns a L x dk(1,1) x nChains Vector
		%
		% by default:
		% {'bath1correlators'}            : projects onto diabatic states
		% {'bath1correlators','adiabatic'}: projects onto adiabatic states
		
		NC = para.nChains;
		out = zeros(para.L, para.dk(1,1), NC);
			
		if para.useStarMPS
			for mc = 1:NC
				% extract MPS for each chain and calculate <anam> separately
				cL = para.chain{mc}.L;
				mpsC  = [mps{1}, cellfun(@(x) x{mc},mps(2:cL),'UniformOutput',false)];
				VmatC = [Vmat{1}, cellfun(@(x) x{mc},Vmat(2:cL),'UniformOutput',false)];
				mpsC{1} = permute(mpsC{1}, [1:mc,mc+2:NC+1, mc+1, NC+2]);
				mpsC{1} = reshape(mpsC{1}, [],para.D(mc,1),para.dk(1,1));

				paraC = para;
				paraC.nChains = 1; paraC.L = para.chain{mc}.L;
				paraC.dk = para.dk(mc,:); paraC.dk(1) = para.dk(1,1);	% need to preserve system dk
				paraC.shift = para.shift(mc,:);
				% calculate <anam>
				if length(type) == 1 || strcmpi(type{2},'diabatic')
					% diabatic projections (H0 eigenstates)
					for kk = 1:paraC.dk(1,1)
						stateProj = zeros(paraC.dk(1,1)); stateProj(kk,kk) = 1;
						bondProj  = [];
						out(1:paraC.L,kk,mc) = calBath1SiteCorrelators_MC(mpsC,VmatC,paraC,stateProj,bondProj);
					end
				elseif strcmpi(type{2},'adiabatic')
					% adiabatic states
					for kk = 1:paraC.dk(1,1)
						stateProj = [];
						bondProj  = zeros(para.D(mc,1)); bondProj(kk,kk) = 1;
						out(1:paraC.L,kk,mc) = calBath1SiteCorrelators_MC(mpsC,VmatC,paraC,stateProj,bondProj);
					end
				end
			end
		else
			if length(type) == 1 || strcmpi(type{2},'diabatic')
				% diabatic projections (H0 eigenstates)
				for kk = 1:para.dk(1,1)
					stateProj = zeros(para.dk(1,1)); stateProj(kk,kk) = 1;
					bondProj  = [];
					out(1:para.L,kk,:) = calBath1SiteCorrelators_MC(mps,Vmat,para,stateProj,bondProj);
				end
			elseif strcmpi(type{2},'adiabatic')
				% adiabatic states
				for kk = 1:para.dk(1,1)
					stateProj = [];
					bondProj  = zeros(para.dk(1,1)); bondProj(kk,kk) = 1;
					out(1:para.L,kk,:) = calBath1SiteCorrelators_MC(mps,Vmat,para,stateProj,bondProj);
				end
			end
% 			out(:,1,:) = calBath1SiteCorrelators_MC(mps,Vmat,para,[1,0;0,0]);	% spin up
% 			out(:,2,:) = calBath1SiteCorrelators_MC(mps,Vmat,para,[0,0;0,1]);	% spin down
		end

	case 'starpolaron'
		%% does mapping chain -> star
		% uses bath1correlators
		% x stands for continuous variable (momentum k)
		% Does state projection!
		% Works with Single and Multi Chain Vmat / Vtens
		%
		%
		if para.foldedChain == 1, error('FoldedChain not yet implemented for polaron'); end
		NC = para.nChains;
		% Get correlators
		if length(type) == 1 || strcmpi(type{2},'diabatic')
			An = real(getObservable({'bath1correlators'},mps,Vmat,para));		% L x states x nChains
		elseif strcmpi(type{2},'adiabatic')
			An = real(getObservable({'bath1correlators','adiabatic'},mps,Vmat,para));		% L x states x nChains
		end
% 		AnUp   = real(calBath1SiteCorrelators_MC(mps,Vmat,para,1));			% L x nChains
% 		AnDown = real(calBath1SiteCorrelators_MC(mps,Vmat,para,-1));		% L x nChains

		[x,~,~] = getStarMapping(para,1);
		out = zeros(size(An,2)+1, length(x), NC);
		for mc = 1:NC % for each chain:
			[x,h,pxn] = getStarMapping(para,mc);
			if strcmp(para.chain{mc}.mapping, 'LanczosTriDiag')
				result = 2.* pxn * An(:,:,mc);	% pxn == U map -> no need for hsquared
				out(:,1:length(x),mc) = [x,result(2:end,:)].';
			else
% 				out(:,:,mc) = [x, 2.*h.*(pxn*An(:,1,mc)),2.*h.*(pxn*An(:,2,mc))]';
				out(:,:,mc) = [x, 2.*(h*ones(1,size(An,2))).*(pxn*An(:,:,mc))]';
			end

		end

    case 'energy'
        % Calculates the entire energy of the chain
		% only usable for focus on spin site!
		% needs type{2} = op
		if length(type) == 2
            out = calEnergy(mps,Vmat,para,type{2});
		else
            error('Need op as 2nd argument');
		end

	case 'coherence'
        % only for spin in SBM

	case 'hshi'
        % Calculates the entire energy of the system hamiltonian and first interaction terms
		% needs type{2} = op
		if length(type) == 2
            out = calHsHi(mps,Vmat,para,type{2});
		else
            op.h1term = cell(1,2); % One body terms of the Hamiltonian
			op.h2term = cell(para.M, 2,2);			% default: nearest neighbor interaction

			for s=1:2
				op=genh1h2term_onesite(para,op,s);
			end
			out = calHsHi(mps,Vmat,para,op);
		end
end

end

function spin = calSpin(mps,Vmat,para)
% Calculate the spin expectation value
% Has to be modified if site 1 changes dimension!!
% Modified:
%   FS 24/05/2014:  excluded MLSpinBoson. Use with try-catch block
%   FS 19/08/2015: - Can be used with Multi-Chain Vmat and Vtens code
%
%   TODO: Needs extension to use information about dimension of spin site.
%
if para.useVtens || para.useStarMPS               % only Quick Fix for Vtens code, Could be used as Standard code!
	McOp = cell(1,1,3);        % L x NC x N
	[McOp{1}, McOp{2}, McOp{3}] = spinop(para.spinbase);
	if para.useVtens
		spinVal = real(expectation_allsites_MC(McOp,mps,Vmat,para));
	else
		spinVal = real(expectation_allsites_StarMPS(McOp,mps,Vmat,para));
	end
	spin.sx = spinVal(1);
	spin.sy = spinVal(2);
	spin.sz = spinVal(3);
	return;
end
N=para.L;
assert(N==length(mps) && N==length(Vmat));
assert(~strcmp(para.model,'MLSpinBoson'),'not possible for MLSBM');

ndset=cell(1,N);

for j=1:N
    ndset{1,j}=eye(size(Vmat{j},1));
end
sx=ndset;sy=sx;sz=sy;

%debug:
% sx{1,1}

[sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
sx{para.spinposition(end)}=sigmaX;
sy{para.spinposition(end)}=sigmaY;
sz{para.spinposition(end)}=sigmaZ;
if strcmp(para.model,'2SpinPhononModel')
    sx{para.spinposition}=kron(sigmaZ,eye(2));  %measures excitation of site1
    sz{para.spinposition}=kron(eye(2),sigmaZ);  %measures excitation of site2
    %try to find a good way to measure this!
    sy{para.spinposition}=kron(sigmaY,eye(2));  %measures only sy of site1
end
spin.sx=expectationvalue(sx,mps,Vmat,mps,Vmat);
spin.sy=expectationvalue(sy,mps,Vmat,mps,Vmat);
spin.sz=expectationvalue(sz,mps,Vmat,mps,Vmat);
spin.sx=real(spin.sx);
spin.sy=real(spin.sy);
spin.sz=real(spin.sz);

end

function nx = calBosonOcc(mps,Vmat,para,varargin)
% Calculate the boson occupation on x chain
% The operator on the spin site is set to zero
% left or right chain occupation via varargin{1} = 1 or 2
%
% Modified:
%	FS 23/01/2014:	- Introduced '~' to ignore unused returned values
%					- support for folded Chain models
%   FS 10/03/2014:  - Using correlator_allsites which is more general
%   FS 04/06/2014:  - updated for spinposition array.
%	FS 17/07/2015:	- varargin = {1,2} to calculate nx/nz for folded chain

n_op = cell(1,para.L);

for j=1:para.L
    if prod(j~=para.spinposition)
        if para.foldedChain == 0
            % 1-chain SBM
            [~,~,n_op{j}] = bosonop(para.dk(j),para.shift(j),para.parity);
        %Modification for 2chain model!! Not perfect or right yet!
        elseif para.foldedChain == 1
            % only kron(n,1) chain occupation calculated.
            if para.parity == 'n'
                [~,~,n] = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
                idm = eye(size(n));
				if varargin{1} == 1
					nr = kron(n,idm);
				elseif varargin{1} == 2
					nr = kron(idm,n);
				end
            else
                [bp,~,~] = bosonop(para.dk(j),para.shift(j),para.parity);   % Why without sqrt??
                [~,~,nr,~,~,~]=paritykron(bp,para.bosonparity);
            end
            n_op{j} = nr;
        else
            disp('Not Implemented: getObservable, calBosonOcc, foldedchain >1');
        end
    else
        n_op{j}=zeros(para.dk(j));      % don't measure spin.
    end
end

%
nx = expectation_allsites(n_op,mps,Vmat);
nx = real(nx);			% imag(nx) = eps -> neglect

end

function n = calBosonOcc_MC(mps,Vtens,para)
%% calculate the boson occupation for Single / Multi-Chain + Vmat / V-tensor-network
% Zero operator for spin site
% calculate for all chains simultaneously
%
% Created 18/08/2015 by FS
assert(para.foldedChain == 0, 'not for folded chains');
L = para.L; NC = para.nChains;

% create Operator for expectationvalue: need nc^2 each (i,:,j) is one operator [] x [] x n x [] x []
McOp = cell(L, NC, NC);
for mc = 1:NC
	for j = 1:para.chain{mc}.L
		if j ~= para.spinposition                     % works with array
			[~,~,McOp{j,mc,mc}] = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
		else
			McOp{j,mc,mc} = zeros(para.dk(1,j));      % don't measure spin.
		end
	end
end

% if para.useStarMPS
% 	n = real(expectation_allsites_StarMPS(McOp,mps,Vtens,para));
% else
	n = real(expectation_allsites_MC(McOp,mps,Vtens,para));  % (NC x L)
% end

end

function bosonshift = calBosonShift(mps,Vmat,para)
% Calculate boson shift x, x^2, var(x)
% The operator on the spin site is set to zero
% Modified:
%	FS 22/01/2014:	- changed to using para.foldedChain.
%   FS 04/06/2014:  - updated for spinposition array
%

L = para.L;
x_opx  = cell(1,L);
x2_opx = cell(1,L);
for j = 1:L
    if j ~= para.spinposition
        if para.foldedChain == 1
            % Constructs Supersite Operators
            % Only measures kron(x,1) chain part
            [bp,~,n] = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
            idm = eye(size(n));
            bpx = kron(bp,idm); bmx = bpx';         %nx = kron(n,idm); unused
            x_opx{j}  = sqrt(2)/2.*(bpx+bmx);
            x2_opx{j} = x_opx{j}*x_opx{j};
        elseif para.foldedChain == 0
            [bp,bm,~] = bosonop(para.dk(1,j),para.shift(1,j),para.parity);
            x_opx{j}  = sqrt(2)/2.*(bp+bm);
            x2_opx{j} = x_opx{j}*x_opx{j};
        end
    else
        x_opx{j}  = zeros(para.dk(j));
        x2_opx{j} = zeros(para.dk(j));
    end
end
bosonshift.x        = expectation_allsites(x_opx,mps,Vmat);
bosonshift.xsquare  = expectation_allsites(x2_opx,mps,Vmat);
bosonshift.xvariant = sqrt(bosonshift.xsquare-bosonshift.x.^2);
bosonshift.xerror   = mean(abs(para.shift-bosonshift.x));
end

function shift = cal1siteShift(mps,Vmat,para)
%% Calculates shift for a single site, para.sitej!
% for use of shifting procedure in optimizesite.m
% mps and Vmat are single-site matrices
% mps has to be focused

NC = para.nChains; j = para.sitej;
xOp_OBB  = cell(1,NC);			% one x-OBB operator per chain
shift    = zeros(1,NC);			% one value per chain
% Generate operators and contract with OBB
if para.foldedChain == 1
	switch para.model
		case 'SpinDoubleBoson' % calculates shift only for x chain!
			bp = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
			if para.parity == 'n'
				idm = eye(size(bp));
				bpx = kron(bp,idm); bmx = bpx';
			else
				[bpx,bmx,~,~,~,~]=paritykron(bp,para.bosonparity);
			end
			x = sqrt(2)/2*(bpx+bmx);
		case '2SpinPhononModel'
			bp = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
			if para.parity=='n'
				idm=eye(size(bp));
				bpr=kron(bp,idm);   bmr=bpr';	% right chain
% 				bpl=kron(idm,bp);   bml=bpl';	% left chain
			else
				[bpr,bmr,~,~,~,~]=paritykron(bp,para.bosonparity);
			end
			x = sqrt(2)/2*(bpr+bmr);                                % why only evaluate it for right part?
	end
	if para.useVmat
		xOp_OBB{1} = Vmat' * x * Vmat;
	else
		xOp_OBB{1} = x;
	end
else
	for mc = 1:NC
		x     = cell(1,NC);
		[bp,bm,~]   = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
		x{mc}       = sqrt(2)/2*(bp+bm);
		xOp_OBB{mc} = contractMultiChainOBB(Vmat, x, para);						% _(n~',n~)
	end
end
temp = contracttensors(conj(mps),3,[1,2],mps,3,[1,2]);							% _(n~',n~)
for mc = 1:NC
	shift(mc) = real(sum(sum(temp.*xOp_OBB{mc})));								% fast tr(A * B.')
end

end

function reducedDensity = calRDM(mps,Vmat,para,k)
% calculates the reduced density matrix of any single site, mps{k} for size(mps{k})=a_{k-1} x a_k x n_k
%
%   created 24/05/2014 by Florian Schroeder @ University of Cambridge
%
%
% copied from prepare.m:
% does l -> r sweep to create state in local picture of k
para.sweepto = 'r';
if length(k) == 1		% single-site RDM
	for i = 1:k-1
		if para.useVmat==1
			[Vmat{i},V] = prepare_onesiteVmat(Vmat{i},para);			% Vmat = U * S * V' ; Vmat := U; V:= S * V'
			mps{i} = contracttensors(mps{i},3,3,V,2,2);                 % = Ai_{l,r,n} * V'_{p,n}; This contraction is defined differently to the paper.
		end
		[mps{i}, U] = prepare_onesite(mps{i}, para,i);             % SVD(Ai_(l,r,n)) = Ai_(l,m,n) * U_(m,r)
		mps{i+1} = contracttensors(U,2,2,mps{i+1},3,1);                 % U_(m,l) * A(i+1)_(l,r,n)
		para=gennonzeroindex(mps,Vmat,para,i);                          % only if parity not 'n'
		para=gennonzeroindex(mps,Vmat,para,i+1);                        % only if parity not 'n'
	end

	% now in form: Al{1}...Al{k-1} M{k} Ar{k+1}...Ar{L}
	%   with Al = left-normalized, Ar: right-normalized.

	n = ndims(mps{k});
	reducedDensity = contracttensors(mps{k},n,1:n-1,conj(mps{k}),n,1:n-1);		% contract rD_nm = Mk_abn Mk*_abm

	% reducedDensity = Vmat{k} * (reducedDensity * Vmat{k}');
	reducedDensity = contracttensors(reducedDensity,2,2,conj(Vmat{k}),2,2);     % contract rD_nj = rD_nm Vmat*_jm
	reducedDensity = contracttensors(Vmat{k},2,2,reducedDensity,2,1);			% contract rD_ij = Vmat_in rd_nj
elseif length(k) == 2 && k(2) == k(1) + 1
	% 2-site RDM of nearest neighbours!
	% move to k(1)
	for i = 1:k(1)-1
		if para.useVmat==1
			[Vmat{i},V] = prepare_onesiteVmat(Vmat{i},para);			% Vmat = U * S * V' ; Vmat := U; V:= S * V'
			mps{i} = contracttensors(mps{i},3,3,V,2,2);                 % = Ai_{l,r,n} * V'_{p,n}; This contraction is defined differently to the paper.
		end
		[mps{i}, U] = prepare_onesite(mps{i}, para,i);             % SVD(Ai_(l,r,n)) = Ai_(l,m,n) * U_(m,r)
		mps{i+1} = contracttensors(U,2,2,mps{i+1},3,1);                 % U_(m,l) * A(i+1)_(l,r,n)
		para=gennonzeroindex(mps,Vmat,para,i);                          % only if parity not 'n'
		para=gennonzeroindex(mps,Vmat,para,i+1);                        % only if parity not 'n'
	end
	% use ncon, since easy and not performance-critical!
	reducedDensity = ncon({mps{k(1)}, mps{k(2)}, conj(mps{k(1)}), conj(mps{k(2)}), Vmat{k(1)}, Vmat{k(2)}, conj(Vmat{k(1)}), conj(Vmat{k(2)})},...
						  {[1,2,3],   [2,4,5],   [1,6,7],         [6,4,8],         [-2,3],     [-1,5],     [-4,7],           [-3,8]});
end
end

function participation = calParticipation(rdm)
% rdm: a reduced density matrix of a single site
%   e.g. rdm = calRDM(mps,Vmat,para,k)
%
% takes a rdm and calculates the participation ratio

participation = 1/sum(diag(rdm).^2);

end

function tunnelE = calTunnelingEnergy(mps,Vmat,para,op)
%% calculates the tunneling energy in MLSBM of the PPC system
%   <Psi|H0-diag(H0)|Psi>
%
%   could be applied to any system.
%   Interacting System Hamiltonian: HI


HI = cell(1);
if ~exist('op','var') || ~isfield(op,'h1term')
    disp('guessing H0 from para for MLSBM');
    op.h1term{1,1} = Hamiltonian_PPC(para);
end
HI{1} = op.h1term{1,1}-diag(diag(op.h1term{1,1}));
tunnelE = expectationvalue(HI,mps,Vmat,mps,Vmat);

end

function An = calBath1SiteCorrelators(mps,Vmat,para,spinProj)  % DEPRECATED. Use calBath1SiteCorrelators_MC instead!
% calculates Re<(1+-sz)/4 * a_n^+>.
% output is matrix containing all values.
% Only calculate Re, since Im is not needed! -> saves space!
% custom made solution for biggest speedup, exact solution!
% spinProj = +-1 to select up/down projection.

An = zeros(para.L,1);						% initialize results array
bp = cell(1,para.L);						% containing all necessary operators
%% generate all operators:
for j=1:para.L
    if prod(j~=para.spinposition)
        if para.foldedChain == 1
            %% not supported yet
			error('This feature is not yet supported');
        elseif para.foldedChain == 0
            [bp{1,j},~,~] = bosonop(para.dk(j),para.shift(j),para.parity);
        end
	else
		[~,~,sz]=spinop(para.spinbase);
        bp{1,j} = (eye(2)+spinProj*sz)/4;			% additional 1/2 for displacement
    end
end

%% compute all partial contractions
% such that (1+sz/4) gets updated to the right in row 1;
% Trace(contraction for result) = Am(j)
bpContract = cell(1,para.L);
for j = 1:para.L
% 	fprintf('%g-',j);
	if j ~= 1
		An(j,1) = trace(updateCleft(bpContract{1,j-1},mps{j},Vmat{j},bp{1,j},mps{j},Vmat{j}));	% this is the result!
	end

	if j ~=1
		bpContract{1,j} = updateCleft(bpContract{1,j-1},mps{j},Vmat{j},[],mps{j},Vmat{j});
	else % j == 1
		bpContract{1,1} = updateCleft([],mps{j},Vmat{j},bp{1,j},mps{j},Vmat{j});                % for spin-projection
	end
end
% fprintf('\n');

end

function An = calBath1SiteCorrelators_MC(mps,Vmat,para,stateProj,bondProj)
% calculates Re<(1+-sz)/4 * a_n^+>.
% output is matrix containing all values for all chains (L x nChains)
% Only calculate Re, since Im is not needed! -> saves space!
% custom made solution for biggest speedup, exact solution!
% 	stateProj: equals state-projection operator, normalised!
%	bondProj:  selecting nth-dominant state of bond
%				[]: no selection, take entire bond
%				n : select nth-state
% supports Multi-Chain with Vmat and Vtens, Works with single chain!
%  does NOT support StarMPS input!

assert(para.foldedChain == 0, 'Please use single-chain code for folded Chains');

An     = zeros(para.L, para.nChains);									% initialize results array
bp_OBB = cell( para.L, para.nChains);
%% generate all operators & contract with OBB Vmat / Vtens:
for j = 1:para.L
	para.sitej = j;
	for mc = 1:para.nChains
		bp = cell(1, para.nChains);						            % containing single-site, single-chain operator
		if j ~= para.spinposition
			bp{mc}        = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
			bp_OBB{j, mc} = contractMultiChainOBB(Vmat{j}, bp, para);
		else	% Spin & Multi-Level System
			if ~isempty(stateProj)
				bp_OBB{j,mc} = Vmat{j}' * (stateProj./2) * Vmat{j};				% Spin-site Vmat, 1/2 since f = (a + a^+)/2
			else
				bp_OBB{j,mc} = (Vmat{j}' * Vmat{j})/2;							% = eye(d_opt)/2; remains unused, since using bondProj below!
			end
		end
	end
end

if isempty(stateProj) && ~isempty(bondProj)
	% do SVD from site 1 to 2 in order to get projections onto the dominant system states!
	%	otherwise have projection onto dominant Bath states (from previous sweep to left)
	A = mps{1};
	[D1, D2, d] = size(A);				% d is site dimension
	A = permute(A, [3, 1, 2]);
    A = reshape(A, [d * D1, D2]);		% reshapes to (a1 d),a2
    [B, S, U] = svd2(A);				% Could also use QR decomposition if nargin !=5
	
	DB = size(S, 1);					% new a2 dimension of B should == d
	B = reshape(B, [d, D1, DB]);
	mps{1} = permute(B, [2, 3, 1]);		% not used anymore!
	% do not contract S into mps{2}, to avoid the displacement to be contaminated with the state amplitude!
	mps{2} = contracttensors(U,2,2, mps{2},3,1);
end

%% compute all partial contractions
% such that (1+sz/4) gets updated to the right in row 1;
% Trace(contraction for result) = Am(j)
bpContract = cell(para.L, 1);										% same for all chains! since only 1 common MPS backbone
for j = 1:para.L
	for mc = 1:para.nChains
		if j ~= 1
			An(j,mc)    = trace(updateCleft(bpContract{j-1,1},mps{j},[],bp_OBB{j,mc},mps{j},[]));	% this is the result!
		end
	end

	if j ~=1
		bpContract{j,1} = updateCleft(bpContract{j-1,1},mps{j},[],[],mps{j},[]);				% Vmat is normalised -> leave out!
	else % j == 1
		if isempty(bondProj)
			bpContract{1,1} = updateCleft([],mps{j},[],bp_OBB{j,1},mps{j},[]);					% for state-projection; = (r',r)
		else
			bpContract{1,1} = bondProj/2;														% 1/2 since f = (a + a^+)/2
		end
	end
end

end

function AmAn = calBath2SiteCorrelators(mps,Vmat,para) % DEPRECATED. Use calBath2SiteCorrelators_MC instead!
% calculates <a_m^+ a_n> for any combination, where n>m.
% output is matrix containing all values.
% Only calculate Re, since Im is not needed! -> saves space!
% custom made solution for biggest speedup, exact solution!

AmAn = zeros(para.L);						% initialize results array
bpbm = cell(2,para.L);						% containing all necessary operators
%% generate all operators:
for j=1:para.L
    if prod(j~=para.spinposition)
        if para.foldedChain == 1
            %% not supported yet
			error('This feature is not yet supported');
		elseif para.nChains > 1
			error('This feature is not yet supported');
        elseif para.foldedChain == 0
            [bp,bm,~] = bosonop(para.dk(j),para.shift(j),para.parity);
            bpbm{1,j} = bp;
			bpbm{2,j} = bm;
        end
    else
        bpbm{1,j} = zeros(para.dk(j));			% zero operators for spin sites
        bpbm{2,j} = zeros(para.dk(j));
    end
end

%% compute all partial contractions
% such that <m n> is meeting in the middle -> reduce overhead!
bpContract = cell(para.L);
for j = 1:para.L-1
	fprintf('%g-',j);
	for i = 1:j+2
		if (i < j) && (j <= floor( (para.L+i)/2))
			bpContract{i,j} = updateCleft(bpContract{i,j-1},mps{j},Vmat{j},[],mps{j},Vmat{j});
		elseif i == j % contruct <leftDM * a+
			if j ~=1
				bpContract{i,i} = updateCleft(bpContract{i,j-1},mps{j},Vmat{j},bpbm{1,j},mps{j},Vmat{j});
			else
				bpContract{i,i} = updateCleft(       []        ,mps{j},Vmat{j},bpbm{1,j},mps{j},Vmat{j});
			end
		elseif i == j+1 % construct left-effective DM
			if j ~= 1
				bpContract{i,j} = updateCleft(bpContract{i-1,j-1},mps{j},Vmat{j},[],mps{j},Vmat{j});
			else
				bpContract{i,j} = updateCleft(        []         ,mps{j},Vmat{j},[],mps{j},Vmat{j});
			end
		elseif i == j+2 % a^+ a of same site
			if j ~= 1
				bpContract{i,j} = trace(updateCleft(bpContract{i-2,j-1},mps{j},Vmat{j},bpbm{1,j}*bpbm{2,j},mps{j},Vmat{j}));
			else
				bpContract{i,j} = trace(updateCleft(        []         ,mps{j},Vmat{j},bpbm{1,j}*bpbm{2,j},mps{j},Vmat{j}));
			end
		end
	end
end
L=para.L;
% The last occupation:
bpContract{L+2,L} = trace(updateCleft(bpContract{L,L-1},mps{L},Vmat{L},bpbm{1,L}*bpbm{2,L},mps{L},Vmat{L}));
fprintf('\n');

bmContract = cell(para.L);
for j = para.L:-1:2
	fprintf('%g-',j);
	for i = para.L:-1:j-1
		if (i > j) && (j > ceil( i/2 ) )
			bmContract{i,j} = updateCright(bmContract{i,j+1},mps{j},Vmat{j},[],mps{j},Vmat{j});
		elseif i == j
			if j ~= para.L
				bmContract{i,i} = updateCright(bmContract{i,j+1},mps{j},Vmat{j},bpbm{2,j},mps{j},Vmat{j});
			else
				bmContract{i,i} = updateCright(        []       ,mps{j},Vmat{j},bpbm{2,j},mps{j},Vmat{j});
			end
		elseif i == j-1
			if j ~= para.L
				bmContract{i,j} = updateCright(bmContract{i+1,j+1},mps{j},Vmat{j},[],mps{j},Vmat{j});
			else
				bmContract{i,j} = updateCright(        []         ,mps{j},Vmat{j},[],mps{j},Vmat{j});
			end
		end
	end
end

for j = 1:para.L
	for i = 1:j
		midPoint = floor((j+i)/2);
		if  i ~= j
			AmAn(i,j) = trace(bpContract{i,midPoint} * bmContract{j,midPoint+1}.');
		else
			AmAn(i,j) = bpContract{i+2,j};
		end
	end
end

end

function AmAn = calBath2SiteCorrelators_MC(mps,Vmat,para)
% calculates <a_m^+ a_n> for any combination, where n>m.
% output is matrix containing all values.
% Only calculate Re, since Im is not needed! -> saves space!
% custom made solution for biggest speedup, exact solution!
%
% Only reduction for Multi-Chain: effective left density matrix
assert(para.foldedChain == 0, 'FoldedChain not yet supported');
AmAn = zeros(para.L,para.L,para.nChains);	% initialize results array
bpbm = cell(3,para.L,para.nChains);			% containing all necessary operators, contracted with OBB
leftDM  = cell(para.L,1);					% left-effective Density Matrix for each site; right-effective should always be 1
rightDM = cell(para.L,1);
%% generate all operators and contract OBB
%  also produce left-effective DM
for j = 1:para.L
	para.sitej = j;
	for mc = 1:para.nChains
		bp = cell(1,para.nChains); bm = cell(1,para.nChains); n = cell(1,para.nChains);
		if j ~= para.spinposition
			[bp{mc}, bm{mc}, n{mc}] = bosonop(para.dk(mc,j),para.shift(mc,j),para.parity);
			bpbm{1,j,mc} = contractMultiChainOBB(Vmat{j}, bp, para);
			bpbm{2,j,mc} = contractMultiChainOBB(Vmat{j}, bm, para);
			bpbm{3,j,mc} = contractMultiChainOBB(Vmat{j},  n, para);
		else
			bpbm{1,j,mc} = zeros(para.dk(1,j));			% zero operators for spin sites
			bpbm{2,j,mc} = zeros(para.dk(1,j));			% dk = d_opt
			bpbm{3,j,mc} = zeros(para.dk(1,j));
		end
	end
	if j ~= 1
		leftDM{j} = updateCleft(leftDM{j-1}, mps{j}, [], [], mps{j}, []);
	else
		leftDM{j} = updateCleft(    []     , mps{j}, [], [], mps{j}, []);
	end
end

%% compute all partial contractions
% such that <m n> is meeting in the middle -> reduce overhead!
for mc = 1:para.nChains
	bpContract = cell(para.L);
	for j = 1:para.L-1
		if j == 1 || mod(j,10) == 0
			fprintf('%g-',j);
		end
		for i = 1:j+2
			if (i < j) && (j <= floor( (para.L+i)/2))
				% take a^+ to next effective basis
				bpContract{i,j} = updateCleft(bpContract{i,j-1},mps{j},[],[],mps{j},[]);
			elseif i == j % contruct <leftDM * a+
				if j ~=1
					bpContract{i,i} = updateCleft(leftDM{j-1},mps{j},[],bpbm{1,j,mc},mps{j},[]);
				else
					bpContract{i,i} = updateCleft(    []     ,mps{j},[],bpbm{1,j,mc},mps{j},[]);
				end
			elseif i == j+1 % construct left-effective DM
				% do nothing, leave empty
			elseif i == j+2 % a^+ a of same site
				if j ~= 1
					bpContract{i,j} = trace(updateCleft(leftDM{j-1},mps{j},[],bpbm{3,j,mc},mps{j},[]));
				else
					bpContract{i,j} = trace(updateCleft(     []    ,mps{j},[],bpbm{3,j,mc},mps{j},[]));
				end
			end
		end
	end
	L = para.L;
	% The last occupation:
	bpContract{L+2,L} = trace(updateCleft(leftDM{L-1},mps{L},[],bpbm{3,L,mc},mps{L},[]));
	fprintf('\n');

	bmContract = cell(para.L);
	for j = para.L:-1:2
		if j == 1 || mod(j,10) == 0
			fprintf('%g-',j);
		end
		for i = para.L:-1:j
			if (i > j) && (j > ceil( i/2 ) )
				bmContract{i,j} = updateCright(bmContract{i,j+1},mps{j},[],[],mps{j},[]);
			elseif i == j
				if j ~= para.L
					bmContract{i,i} = updateCright(rightDM{j+1},mps{j},[],bpbm{2,j,mc},mps{j},[]);
				else
					bmContract{i,i} = updateCright(     []     ,mps{j},[],bpbm{2,j,mc},mps{j},[]);
				end
			end
		end
		if mc == 1 % construct right effective DM only once
			if j ~= para.L
				rightDM{j} = updateCright(rightDM{j+1},mps{j},[],[],mps{j},[]);
			else
				rightDM{j} = updateCright(     []     ,mps{j},[],[],mps{j},[]);
			end
		end
	end

	for j = 1:para.L
		for i = 1:j
			midPoint = floor((j+i)/2);
			if  i ~= j
				AmAn(i,j,mc) = sum(sum(bpContract{i,midPoint} .* bmContract{j,midPoint+1}));		% faster tr(A * B.')
			else
				AmAn(i,j,mc) = bpContract{i+2,j};
			end
		end
	end
end
end

function E = calEnergy(mps,Vmat,para,op)
	%%
	sitej=1;
	if para.useStarMPS
		E = 0;			% nut supported for now!
		return;
	end
% 	[~,BondDimRight] = size(op.Hright);
% 	[~,BondDimLeft]  = size(op.Hleft);
% 	op = h1h2toOBB(Vmat{sitej},para,op);		% DEPRECATED
	op  = H_Eff([],Vmat{sitej},'A',op,para);
% 	[~,OBBDim]	 = size(op.h1j);
	[BondDimLeft, BondDimRight, OBBDim] = size(mps{sitej});
	M = size(op.h2j,1);
	A = reshape(mps{sitej},[numel(mps{sitej}),1]);
	E = A'*HmultA(A, op, BondDimLeft, BondDimRight, OBBDim, M,para.parity,[]);
	E = real(E);		% imag(E) = eps -> negligible!
end

function n = expectation_allsites_MC(McOp,mps,Vmat,para)
%% Multi Chain Expectation values for Vmat and Vtens
% copied from expectation_allsites since faster to write
%
% McOp : L x NC x N
% NC = para.nChains;
if para.useStarMPS		% shortcut for StarMPS
	n = expectation_allsites_StarMPS(McOp,mps,Vmat,para);
	return;
end
[L,~,N] = size(McOp);               % N = number of "observable chains"
% assert(L == length(mps));			% would also work for N < dim(mps)

% Contract McOp with OBB of Vmat or Vtens
McOp_OBB = cell(N,L);
for j = 1:L
	para.sitej = j;
	for n = 1:N
		if j ~= para.spinposition
			McOp_OBB{n,j} = contractMultiChainOBB(Vmat{j}, McOp(j,:,n), para);
		else
			Ind = ~cellfun('isempty', McOp(j,:,n));
			McOp_OBB{n,j} = Vmat{j}' * McOp{j,Ind,n} * Vmat{j};
		end
	end
end

% call with empty Vmat since OBB is already contracted
n = expectation_allsites(McOp_OBB, mps, cell(1,para.L));	% N x L

end

function n = expectation_allsites_StarMPS(McOp,mps,Vmat,para)
%% Multi Chain Expectation values for StarMPS code
% copied from expectation_allsites since faster to write
%
% McOp : L x NC x N
NC = para.nChains;
[L,M,N] = size(McOp);               % N = number of "observable chains";
% assert(L == length(mps));			% would also work for N < dim(mps)

n = zeros(N,L);
% Here: iterate over chains and determine whether McOp is empty
for mc = 1:M
	tempOp = reshape(McOp(:,mc,:),L,N);				% L x N
	empty  = sum(cellfun('isempty',tempOp),1);		% 0 where Op, L where empty
	if all(empty == L), continue; end
	% else: extract mpsChain, mpsVmat including system site
	cL = para.chain{mc}.L;
	mpsChain  = [mps(1), cellfun(@(x) x{mc},mps(2:cL),'UniformOutput',false)];
	VmatChain = [Vmat(1), cellfun(@(x) x{mc},Vmat(2:cL),'UniformOutput',false)];
	mpsChain{1} = permute(mpsChain{1}, [1:mc,mc+2:NC+1, mc+1, NC+2]);
	mpsChain{1} = reshape(mpsChain{1}, [],para.D(mc,1),para.dk(1,1));
	for ii = find(empty ~= L)
		n(ii,1:min(L,cL)) = expectation_allsites(tempOp(1:min(L,cL),ii).', mpsChain, VmatChain);
	end
end

end

function [x,hsquared,pxn] = getStarMapping(para,mc)
	%% finds mapping chain -> star
		n = 0:(para.L-2);		% length(n) = L-1
		alphaN   = zeros(length(n),1);
		betaN    = zeros(length(n),1);			% this is sqrt(betaN); betaN(1) = beta_0 is meaningless!

		% find bath parameters:
		if strcmp(para.chain{mc}.mapping,'OrthogonalPolynomials')
			step = max(1/para.L/10,10^-3);
			x = (0:step:1)';		% resolution
			hsquared = ones(length(x),1);
			pxn = zeros(length(x),para.L);			% these are orthonormalized polynomials. examples:
					% pxn(:,1) = 0; pxn(:,2) = px(0) = 1/t(1); pxn(:,3) =
					% x-alpha(0)/norm();
			s = para.chain{mc}.s;
			if strcmp(para.chain{mc}.spectralDensity,'Leggett_Hard')
				hsquared = 2.*para.chain{mc}.alpha.*(x.^s);
				alphaN = (1+( (s^2)./((s+2.*n).*(2+s+2.*n)) ))./2;
				betaN  = ((1+n).*(1+s+n))./(s+2+2.*n)./(3+s+2.*n).*sqrt((3+s+2.*n)./(1+s+2.*n));
			elseif strcmp(para.chain{mc}.spectralDensity,'Leggett_Soft')
				x = x.*3;						% need some longer range since spectralDensity is not limited in bandwidth!
				hsquared = 2.*para.chain{mc}.alpha.*(x.^s).*exp(-x);
				alphaN = 2.*n+1+s;
				betaN  = sqrt((n+1).*(n+1+s));
			end
		elseif strcmp(para.chain{mc}.mapping, 'Stieltjes')
			x        = reshape(para.chain{mc}.xi, [],1);
			hsquared = reshape(para.chain{mc}.Gamma.^2,[],1);
			alphaN   = reshape(para.chain{mc}.epsilon, 1,[]);
			betaN    = reshape(para.chain{mc}.t(2:end),1,[]);
			pxn      = zeros(length(x),para.chain{mc}.L);
		elseif strcmp(para.chain{mc}.mapping, 'LanczosTriDiag')
			x        = reshape(para.chain{mc}.xi, [],1);			% frequency Axis
			hsquared = reshape(para.chain{mc}.Gamma.^2,[],1);		% not needed
			pxn      = para.chain{mc}.U;							% not poly, but map U
			return;
		end
		% find OrthPol:
		pxn(:,2) = ones(length(x),1).*1/para.chain{mc}.t(1);		% p(0) = constant normalized pol.
		for i = 1:para.L-2
			% p(i+1) = ((x - a(i))p(i)-b(i)p(i-1))/b(i+1)
			pxn(:,i+2) = ( (x - alphaN(i)).*pxn(:,i+1) - betaN(i).*pxn(:,i) )./betaN(i+1);
		end


end

function E = calHsHi(mps,Vmat,para,op)
%% function E = calHsHi(mps,Vmat,para,op)
% calculates the expectation value of <H_s + H_i(1,2)>
% of only System and first bond.
if ~exist('op','var') || ~isfield(op,'h1term')
	error('op is needed to use calHsHi');
end

% create n_op for expectation_general_StarMPS
% n_op = M+1 x L x NC
%	where H1term = n_op(1,1,1)
M  = para.M;
NC = para.nChains;
L  = para.L;

n_op = cell(M+1,L,NC);

n_op{1,1,1} = op.h1term{1,1};
for m = 1:M
	idx = squeeze(~cellfun('isempty',op.h2term(m,2,2,:)));
	n_op(m+1,1,idx) = op.h2term(m,1,1,1);
	n_op(m+1,2,idx) = op.h2term(m,2,2,idx);
end

E_i = expectation_general_StarMPS(n_op,mps,Vmat,para);
E = [real(E_i(1));2*real(E_i(2:2:end))]';				% find a better solution for this! But not now..

% E = sum(E_i);

end

function n = expectation_general(n_op, mps, Vmat, para)
%% function n = expectation_general(n_op, mps, Vmat, para)
% computes the expectation value of an arbitrary (general) multi-site operator defined as:
% O = 1 x L cell
%		where site-operators are stored in respective cells. Other cells: []
%	e.g:  < O_i O_j O_k >
%		O = {[],...,[],O_i,[],...,[],O_j,[],...,[],O_k,[],...}
%
% if multiple operators are desired on the same chain
% n_Op = M x L cell
%       where n_Op{1,:} = O{1}
%             n_Op{2,:} = O{2}  ...
%
% only for single-chain MPS
% output contains no site-dependence, such as expectation_allsites(). This has to be dealt with inside the calling function!
%
% written by FS 09/11/2015

[M,~] = size(n_op);						% M x L
empty = cellfun('isempty',n_op);

tempC = cell(M,1);						% will contain temporary right-effective operators.
lastnz = find(sum(~empty,1),1,'last');	% last non-zero

for ii = lastnz:-1:1
	% sweep from back
	for m = 1:M
		if ii == 1 && size(mps{1},1) > size(mps{1},2)
			n(m) = trace(updateCleft([],mps{ii},Vmat{ii},n_op{m,ii},mps{ii},Vmat{ii})*tempC{m}.');		% contract left first, since memory can blow up!
		else
			tempC{m} = updateCright(tempC{m},mps{ii},Vmat{ii},n_op{m,ii},mps{ii},Vmat{ii});
			if ii == 1
				n(m) = trace(tempC{m});
			end
		end
	end
end

end

function n = expectation_general_StarMPS(n_op, mps, Vmat, para)
%% function n = expectation_general(n_op, mps, Vmat, para)
% computes the expectation value of an arbitrary (general) multi-site operator defined as:
% O = 1 x L cell
%		where site-operators are stored in respective cells. Other cells: []
%	e.g:  < O_i O_j O_k >
%		O = {[],...,[],O_i,[],...,[],O_j,[],...,[],O_k,[],...}
%
% if multiple operators are desired on any chain mc
% n_Op = M x L x NC cell
%       where n_Op{1,:,mc} = O{1}
%             n_Op{2,:,mc} = O{2}  ...
%
%
% only for single-chain MPS
% output contains no site-dependence, such as expectation_allsites(). This has to be dealt with inside the calling function!
%
% written by FS 09/11/2015

[M,~,NC] = size(n_op);					% M x Lmax x NC
empty = cellfun('isempty',n_op);

n = zeros(M,1);

for mc = 1:NC
	% find which m can be computed
	allnz = find(sum(~empty(:,:,mc),2));	% all non-zero for mc

	% reshape MPS into single-chain
	cL = para.chain{mc}.L;
	mpsChain  = [mps{1}, cellfun(@(x) x{mc},mps(2:cL),'UniformOutput',false)];
	VmatChain = [Vmat{1}, cellfun(@(x) x{mc},Vmat(2:cL),'UniformOutput',false)];
	mpsChain{1} = permute(mpsChain{1}, [1:mc,mc+2:NC+1, mc+1, NC+2]);
	mpsChain{1} = reshape(mpsChain{1}, [],para.D(mc,1),para.dk(1,1));

	% calculate the observables for that chain:
	n(allnz) = expectation_general(n_op(allnz,:,mc),mpsChain,VmatChain,para);
end


end