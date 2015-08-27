function op = genh1h2term_onesite(para,op,s)
%Define the hamiltonian
% op.h2term{i,j,k,l}:
%	l is chain number (only used in SpinBosonMC)
%   k is site number
%   i is number of term in sum to address:  t1*a*b + t1*b*a
%                       sum:                    1       2
%   j is position nr in each hopping term:  t*a*b
%       t1*a*b:     para.t(1)*op.h2term{1,1,1}*op.h2term{1,2,2}
%       t1*b*a:     para.t(1)*op.h2term{2,1,1}*op.h2term{2,2,2}
%example: op.h2term{1,1,1} is with op.h2term{1,2,2}
% interaction terms don't seem to be normal ordered...
%
% If Model changed, modify also:
%   VMPS_SBM1.m:    para.dk(1) gives dimension of first site;
%   calspin.m:      sx,sy,sz defines spin measure operator. Change if dim(site1) changes!
switch para.model
    case 'SpinBoson'
        %%%%%%%%%%%%%%%%%%%Spin-boson Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch s
            case 1                                                  % first chain pos = all spin sites!
                [sigmaX,~,sigmaZ] = spinop(para.spinbase);       % gives XYZ operators with respect to specific main base
                zm_spin			  = zeros(2);
                op.h1term{1}	  = -para.hx./2.*sigmaX-para.hz./2.*sigmaZ;
                op.h2term{1,1,1}  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{1,2,1} = zm_spin;		% t(1) = sqrt(eta_0/pi)/2
                op.h2term{2,1,1}  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{2,2,1} = zm_spin;
            case para.L                                             % last chain pos: only one coupling?
                [bp,bm,n]		  = bosonop(para.dk(para.L),para.shift(para.L),para.parity);
                zm				  = sparse(size(bp,1),size(bp,1));
                op.h1term{s}	  = para.chain{1}.epsilon(para.L-1).*n;
                op.h2term{1,1,s}  = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,s}  = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                [bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
                op.h1term{s}	  = para.chain{1}.epsilon(s-1).*n;                          % e(1) == w(0)
                op.h2term{1,1,s}  = para.chain{1}.t(s).*bp; op.h2term{1,2,s} = bm;    % t(2) == t(n=0)
                op.h2term{2,1,s}  = para.chain{1}.t(s).*bm; op.h2term{2,2,s} = bp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'SpinBosonTTM'
        %%%%%%%%%%%%%%%%%%% Spin-boson Model for the Transfer Tensor Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch s
			case 1                                                  % first chain pos = ancilla system, needs to be maximally entangled with site 2
                zm_spin			  = zeros(2);
                op.h1term{s}	  = zm_spin;
                op.h2term{1,1,s}  = zm_spin; op.h2term{1,2,s} = zm_spin;		% t(1) = sqrt(eta_0/pi)/2
                op.h2term{2,1,s}  = zm_spin; op.h2term{2,2,s} = zm_spin;
            case 2                                                  % first chain pos = all spin sites!
                [sigmaX,~,sigmaZ] = spinop(para.spinbase);       % gives XYZ operators with respect to specific main base
                zm_spin			  = zeros(2);
                op.h1term{s}	  = -para.hx./2.*sigmaX-para.hz./2.*sigmaZ;
                op.h2term{1,1,s}  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{1,2,s} = zm_spin;		% t(1) = sqrt(eta_0/pi)/2
                op.h2term{2,1,s}  = para.chain{1}.t(1).*sigmaZ./2; op.h2term{2,2,s} = zm_spin;
            case para.L                                             % last chain pos: only one coupling?
                [bp,bm,n]		  = bosonop(para.dk(para.L),para.shift(para.L),para.parity);
                zm				  = sparse(size(bp,1),size(bp,1));
                op.h1term{s}	  = para.chain{1}.epsilon(para.L-1).*n;
                op.h2term{1,1,s}  = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,s}  = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                [bp,bm,n]		  = bosonop(para.dk(s),para.shift(s),para.parity);
                op.h1term{s}	  = para.chain{1}.epsilon(s-2).*n;                          % e(1) == w(0)
                op.h2term{1,1,s}  = para.chain{1}.t(s-1).*bp; op.h2term{1,2,s} = bm;    % t(2) == t(n=0)
                op.h2term{2,1,s}  = para.chain{1}.t(s-1).*bm; op.h2term{2,2,s} = bp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'SpinBoson2folded'
        %%%%%%%%%%%%%%%%%%%2-chain Spin Boson Model - One Super-Chain%%%%%%%%%%%%%%%%%%%%%%
		% Created 07/07/15 by F.S.
        switch s
            case 1
                [sigmaX,~,sigmaZ]=spinop(para.spinbase);
                zm_spin=zeros(2);
                op.h1term{1}     = - para.hx./2.*sigmaX - para.hz./2.*sigmaZ;
                op.h2term{1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{1,2,1} = zm_spin;	% X chain
                op.h2term{2,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{2,2,1} = zm_spin;
                op.h2term{3,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{3,2,1} = zm_spin;	% Z chain
                op.h2term{4,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{4,2,1} = zm_spin;
            case para.L
                [bp,~,n] = bosonop(sqrt(para.dk(para.L)),para.shift(para.L),para.parity);  % gives [bp, bm, n]
                if para.parity=='n'
                    idm=eye(size(n));
                    bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
                    bpz=kron(idm,bp);bmz=bpz';nz=kron(idm,n);
				else
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
%                     [bpx,bmx,nx,bpy,bmy,ny]=paritykron(bp,para.bosonparity);
                end
                zm=sparse(size(bpx,1),size(bpx,1));
                op.h1term{para.L}     = para.chain{1}.epsilon(para.L-1).*nx + para.chain{2}.epsilon(para.L-1).*nz;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bmx;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bpx;
                op.h2term{3,1,para.L} = zm; op.h2term{3,2,para.L} = bmz;
                op.h2term{4,1,para.L} = zm; op.h2term{4,2,para.L} = bpz;
            otherwise
                [bp,~,n] = bosonop(sqrt(para.dk(s)),para.shift(s),para.parity);
                if para.parity=='n'
                    idm=eye(size(n));
                    bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
                    bpz=kron(idm,bp);bmz=bpz';nz=kron(idm,n);
				else
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
%                     [bpx,bmx,nx,bpy,bmy,ny]=paritykron(bp,para.bosonparity);
                end
%                 zm=sparse(size(bpx,1),size(bpx,1));
                op.h1term{s}     = para.chain{1}.epsilon(s-1).*nx + para.chain{2}.epsilon(s-1).*nz;
                op.h2term{1,1,s} = para.chain{1}.t(s).*bpx; op.h2term{1,2,s} = bmx;
                op.h2term{2,1,s} = para.chain{1}.t(s).*bmx; op.h2term{2,2,s} = bpx;
                op.h2term{3,1,s} = para.chain{2}.t(s).*bpz; op.h2term{3,2,s} = bmz;
                op.h2term{4,1,s} = para.chain{2}.t(s).*bmz; op.h2term{4,2,s} = bpz;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'SpinDoubleBoson'
        %%%%%%%%%%%%%%%%%%%Spin double boson Model One Chain%%%%%%%%%%%%%%%%%%%%%%
		% From Cheng, DEPRECATED, para.foldedChain = 1!
		% 2 Boson chains folded into supersites
        switch s
            case 1
                [sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
                zm_spin=zeros(2);
                %assert(para.Delta==0);
                op.h1term{1}=-para.hx./2.*sigmaX-para.hy./2.*sigmaY-para.hz./2.*sigmaZ;
                op.h2term{1,1,1} = para.t(1).*sigmaX./2; op.h2term{1,2,1} = zm_spin;
                op.h2term{2,1,1} = para.t(1).*sigmaX./2; op.h2term{2,2,1} = zm_spin;   %to the right
                op.h2term{3,1,1} = para.t(1).*sigmaY./2; op.h2term{3,2,1} = zm_spin;
                op.h2term{4,1,1} = para.t(1).*sigmaY./2; op.h2term{4,2,1} = zm_spin;   %to the left
            case para.L
                [bp,~,n] = bosonop(sqrt(para.dk(para.L)),para.shift(para.L),para.parity);  % gives [bp, bm, n]
                if para.parity=='n'
                    idm=eye(size(n));
                    bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
                    bpy=kron(idm,bp);bmy=bpy';ny=kron(idm,n);
                else
                    [bpx,bmx,nx,bpy,bmy,ny]=paritykron(bp,para.bosonparity);
                end
                zm=sparse(size(bpx,1),size(bpx,1));
                op.h1term{para.L}=para.epsilon(para.L-1).*nx+para.epsilon(para.L-1).*ny;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bmx;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bpx;
                op.h2term{3,1,para.L} = zm; op.h2term{3,2,para.L} = bmy;
                op.h2term{4,1,para.L} = zm; op.h2term{4,2,para.L} = bpy;
            otherwise
                [bp,~,n] = bosonop(sqrt(para.dk(s)),para.shift(s),para.parity);
                if para.parity=='n'
                    idm=eye(size(n));
                    bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
                    bpy=kron(idm,bp);bmy=bpy';ny=kron(idm,n);
                else
                    [bpx,bmx,nx,bpy,bmy,ny]=paritykron(bp,para.bosonparity);
                end
%                 zm=sparse(size(bpx,1),size(bpx,1));
                op.h1term{s}=para.epsilon(s-1).*nx+para.epsilon(s-1).*ny;
                op.h2term{1,1,s} = para.t(s).*bpx; op.h2term{1,2,s} = bmx;
                op.h2term{2,1,s} = para.t(s).*bmx; op.h2term{2,2,s} = bpx;
                op.h2term{3,1,s} = para.t(s).*bpy; op.h2term{3,2,s} = bmy;
                op.h2term{4,1,s} = para.t(s).*bmy; op.h2term{4,2,s} = bpy;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'SpinDoulbeBosonFolded' %%Two chain case. Doesn't work any more.
        %%%%%%%%%%%%%%%%%%Spin Doulbe Boson Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % From Cheng, DEPRECATED
		switch s
            case para.spinposition      %at the middle of the chain?
                if para.parity=='n'
                    sigmaX=[0 1;1 0]; %In the non parity basis
                    sigmaY=[0 -1i;1i 0];
                    sigmaZ=[1 0;0 -1];
                else
                    sigmaX=[-1 0;0 1]; %In the parity basis
                    sigmaZ=[0 -1;-1 0];
                end
                op.h1term{s}=-para.Delta./2.*sigmaX+para.epsilon_spin./2.*sigmaZ;
                op.h2term{1,1,s} = para.tz(1).*sigmaZ./2; op.h2term{1,2,s} = para.ty(1).*sigmaY./2;
                op.h2term{2,1,s} = para.tz(1).*sigmaZ./2; op.h2term{2,2,s} = para.ty(1).*sigmaY./2;
            case 1
                [bp,bm,n] = bosonop(para.dk(1),para.shift(1),para.parity);
                zm=sparse(size(bp,1),size(bp,1));
                op.h1term{1}=para.epsilony(para.Ly-1).*n;
                op.h2term{1,1,1} = bm; op.h2term{1,2,1} = zm;
                op.h2term{2,1,1} = bp; op.h2term{2,2,1} = zm;
            case para.L
                [bp,bm,n] = bosonop(para.dk(para.L),para.shift(para.L),para.parity);
                zm=sparse(size(bp,1),size(bp,1));
                op.h1term{para.L}=para.epsilonz(para.Lz-1).*n;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                if s<para.spinposition
                    ss=para.spinposition-s;
                    [bp,bm,n] = bosonop(para.dk(s),para.shift(s),para.parity);
                    op.h1term{s}=para.epsilony(ss).*n;
                    op.h2term{1,1,s} = bm; op.h2term{1,2,s} = para.ty(ss+1).*bp;
                    op.h2term{2,1,s} = bp; op.h2term{2,2,s} = para.ty(ss+1).*bm;
                else
                    ss=s-para.spinposition;
                    [bp,bm,n] = bosonop(para.dk(s),para.shift(s),para.parity);
                    op.h1term{s} = para.epsilonz(ss).*n;
                    op.h2term{1,1,s} = para.tz(ss+1).*bp; op.h2term{1,2,s} = bm;
                    op.h2term{2,1,s} = para.tz(ss+1).*bm; op.h2term{2,2,s} = bp;
                end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'DissipativeHamonicOscillator'
        %%%%%%%%%%%%%%%%%%%Dissipative Hamonic Oscillator Model%%%%%%%%%%%%%%%%%%%
        switch s
            case 1
                [bp,bm,n] = bosonop(para.dk(1),para.shift(1));
                zm=sparse(size(bp,1),size(bp,1));
                op.h1term{1}=para.Omega.*n+para.epsilon_osc./2.*(bp+bm);
                op.h2term{1,1,1} = para.t(1).*(bp+bm)./2; op.h2term{1,2,1} = zm;
                op.h2term{2,1,1} = para.t(1).*(bp+bm)./2; op.h2term{2,2,1} = zm;
            case para.L
                [bp,bm,n] = bosonop(para.dk(para.L),para.shift(para.L));
                zm=sparse(size(bp,1),size(bp,1));
                op.h1term{para.L}=para.epsilon(para.L-1).*n;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                [bp,bm,n] = bosonop(para.dk(s),para.shift(s));
                op.h1term{s}=para.epsilon(s-1).*n;
                op.h2term{1,1,s} = para.t(s).*bp; op.h2term{1,2,s} = bm;
                op.h2term{2,1,s} = para.t(s).*bm; op.h2term{2,2,s} = bp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case '2SpinPhononModel'
        %%%%%%%%%%%%%%%%%%%2-State Phonon Model (folded)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % copied from SpinDoubleBoson. This is indeed the folded model, as the bosonic sites have kron()
		% possibly not working!
        switch s
            case 1      % if I want to include two exciton sites with chain-ex1-ex2-chain I need to have ex1,ex2 in site 1 combined!
                        % start off with only one site and extend it for later!
                        % possibly wrong: COUPLING in Prior, Chin 10: parallel to energy splitting! so sigmaZ
                [sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
                sigmaPlus = (sigmaX+1i.*sigmaY)./2;
                sigmaMinus = (sigmaX - 1i.*sigmaY)./2;
                zm_spin=zeros(4);
                %assert(para.Delta==0);
                %op.h1term{1}=-para.hx./2.*sigmaX-para.hy./2.*sigmaY-para.hz./2.*sigmaZ;
                op.h1term{1}= para.hx./2.*kron(sigmaZ,eye(2)) + para.hz./2.*kron(eye(2),sigmaZ) + ...
                              para.hy.*(kron(sigmaPlus,sigmaMinus) + kron(sigmaMinus,sigmaPlus));
                %here: hx=epsilon1; hz=epsilon2; hy = J for Prior 2010 with kron(site1,site2);

                op.h2term{1,1,1} = para.t(1)./4.*kron(eye(2)+sigmaZ,eye(2)); op.h2term{1,2,1} = zm_spin;
                op.h2term{2,1,1} = para.t(1)./4.*kron(eye(2)+sigmaZ,eye(2)); op.h2term{2,2,1} = zm_spin; %to the right
                op.h2term{3,1,1} = para.t(1)./4.*kron(eye(2),eye(2)+sigmaZ); op.h2term{3,2,1} = zm_spin;
                op.h2term{4,1,1} = para.t(1)./4.*kron(eye(2),eye(2)+sigmaZ); op.h2term{4,2,1} = zm_spin; %to the left
            case para.L
                [bp,~,n] = bosonop(sqrt(para.dk(para.L)),para.shift(para.L),para.parity);
                if para.parity=='n'     % always use n (might be easier?)
                    idm=eye(size(n));
                    bpr=kron(bp,idm);   bmr=bpr';   nr=kron(n,idm);     % matrices for right    (=x)
                    bpl=kron(idm,bp);   bml=bpl';   nl=kron(idm,n);     % matrices for left     (=y)
                else
                    [bpr,bmr,nr,bpl,bml,nl]=paritykron(bp,para.bosonparity);
                end
                zm=sparse(size(bpr,1),size(bpr,1));
                op.h1term{para.L}=para.epsilon(para.L-1).*nr+para.epsilon(para.L-1).*nl;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bmr;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bpr;
                op.h2term{3,1,para.L} = zm; op.h2term{3,2,para.L} = bml;
                op.h2term{4,1,para.L} = zm; op.h2term{4,2,para.L} = bpl;
            otherwise
                [bp,~,n] = bosonop(sqrt(para.dk(s)),para.shift(s),para.parity);
                if para.parity=='n'
                    idm=eye(size(n));
                    bpr=kron(bp,idm);   bmr=bpr';   nr=kron(n,idm);
                    bpl=kron(idm,bp);   bml=bpl';   nl=kron(idm,n);
                else
                    [bpr,bmr,nr,bpl,bml,nl]=paritykron(bp,para.bosonparity);
                end
                zm=sparse(size(bpr,1),size(bpr,1));
                op.h1term{s}=para.epsilon(s-1).*nr+para.epsilon(s-1).*nl;
                op.h2term{1,1,s} = para.t(s).*bpr; op.h2term{1,2,s} = bmr;
                op.h2term{2,1,s} = para.t(s).*bmr; op.h2term{2,2,s} = bpr;
                op.h2term{3,1,s} = para.t(s).*bpl; op.h2term{3,2,s} = bml;
                op.h2term{4,1,s} = para.t(s).*bml; op.h2term{4,2,s} = bpl;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'MLSpinBoson'
        %%%%%%%%%%%%%%%%%%%Multi-Level Spin-boson Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for use of this model:
        %   para.hx: de
        %   HSI is a matrix defining coupling of system to bath. Is scaled with para.t(1)
        switch s
            case 1                                                  % first chain pos = all spin sites!
                [HS0, HSI] = MLSB_Operators(para);                  % uses local function
                zm_spin=zeros(para.MLSB_Ls);
                op.h1term{1}= HS0;
                op.h2term{1,1,1} = para.t(1).*HSI; op.h2term{1,2,1} = zm_spin;		% t(1) = sqrt(eta_0/pi) only for J_Renger!
                op.h2term{2,1,1} = para.t(1).*conj(HSI'); op.h2term{2,2,1} = zm_spin;		% was sigmaX; conj() to still have hermitian Hamiltonian!
            case para.L
                [bp,bm,n] = bosonop(para.dk(para.L),para.shift(para.L),para.parity);
                zm=sparse(size(bp,1),size(bp,1));
                op.h1term{para.L}=para.epsilon(para.L-1).*n;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                [bp,bm,n] = bosonop(para.dk(s),para.shift(s),para.parity);
                op.h1term{s}=para.epsilon(s-1).*n;
                op.h2term{1,1,s} = para.t(s).*bp; op.h2term{1,2,s} = bm;
                op.h2term{2,1,s} = para.t(s).*bm; op.h2term{2,2,s} = bp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'ImpurityQTN'
        %%%%%%%%%%%%%%%%%%%Impurity coupled to Quantum Telegraph Noise Model%%%%%%%%%%%%%%%%%%%
        % for use of this model:
        %   TODO: this is by far not completed!!! and ready to use..
        %
		switch s
            case 1                                                  % first chain pos = all spin sites!
                [sigmaX,~,sigmaZ]=spinop(para.spinbase);            % gives XYZ operators with respect to specific main base
                zm_spin=zeros(2);
                op.h1term{1}=para.hz./2.*sigmaZ;
                op.h2term{1,1,1} = para.t(1).*sigmaZ./2; op.h2term{1,2,1} = zm_spin;		% was sigmaX
                op.h2term{2,1,1} = para.t(1).*sigmaZ./2; op.h2term{2,2,1} = zm_spin;		% was sigmaX
            case para.L                                             % last chain pos: only one coupling?
                [bp,bm,n] = bosonop(para.dk(para.L),para.shift(para.L),para.parity);
                zm=sparse(size(bp,1),size(bp,1));
                op.h1term{para.L}=para.epsilon(para.L-1).*n;
                op.h2term{1,1,para.L} = zm; op.h2term{1,2,para.L} = bm;
                op.h2term{2,1,para.L} = zm; op.h2term{2,2,para.L} = bp;
            otherwise
                [bp,bm,n] = bosonop(para.dk(s),para.shift(s),para.parity);
                op.h1term{s}=para.epsilon(s-1).*n;
                op.h2term{1,1,s} = para.t(s).*bp; op.h2term{1,2,s} = bm;
                op.h2term{2,1,s} = para.t(s).*bm; op.h2term{2,2,s} = bp;
		end

	case 'SpinBoson2C'
        %%%%%%%%%%%%%%%%%%% Spin-Boson Model - 2-Chain %%%%%%%%%%%%%%%%%%%%%%
		% Not linear, but in multi-chain configuration!
		% Benchmark for multi-chain method! Should yield same as
		%		'SpinBoson2folded'
		% Spin is always in chain 1 for backward compatibility
		%
		% working (23/07/15)
		% Created 15/07/15 by F.S.
        switch s
            case 1
				[sigmaX,~,sigmaZ]  = spinop(para.spinbase);
                zm_spin			   = zeros(2);
                op.h1term{1,1}     = - para.hx./2.*sigmaX - para.hz./2.*sigmaZ;
                op.h2term{1,1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{1,2,1,1} = zm_spin;	% X chain
                op.h2term{2,1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{2,2,1,1} = zm_spin;
                op.h2term{3,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{3,2,1,1} = zm_spin;	% Z chain
                op.h2term{4,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{4,2,1,1} = zm_spin;
            case para.L
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
				[bpz,bmz,nz] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
                if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
                end
                zmx = sparse(size(bpx,1),size(bpx,1));
				zmz = sparse(size(bpz,1),size(bpz,1));
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = para.chain{2}.epsilon(s-1).*nz;
                op.h2term{1,1,s,1} = zmx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = zmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = zmz; op.h2term{3,2,s,2} = bmz;
                op.h2term{4,1,s,2} = zmz; op.h2term{4,2,s,2} = bpz;
            otherwise
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
				[bpz,bmz,nz] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
				if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
				end
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = para.chain{2}.epsilon(s-1).*nz;
                op.h2term{1,1,s,1} = para.chain{1}.t(s).*bpx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = para.chain{1}.t(s).*bmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = para.chain{2}.t(s).*bpz; op.h2term{3,2,s,2} = bmz;
                op.h2term{4,1,s,2} = para.chain{2}.t(s).*bmz; op.h2term{4,2,s,2} = bpz;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'SpinBoson2CT'
        %%%%%%%%%%%%%%%%%%% Spin-Boson Model - 2-Chain Thermal Bath %%%%%%%%%%%%%%%%%%%%%%
		% SBM including ancilla chain for thermal bath states
		% Multi-chain configuration for Vtens
		%
		% Spin is always in chain 1
		%
		% working 27/08/15
		% Created 25/08/15 by F.S.
        switch s
            case 1
				[sigmaX,~,sigmaZ]  = spinop(para.spinbase);
                zm_spin			   = zeros(2); zm = []; %zeros(2);
                op.h1term{1,1}     = - para.hx./2.*sigmaX - para.hz./2.*sigmaZ;
                op.h2term{1,1,1,1} = para.chain{1}.t(1).*sigmaZ./2; op.h2term{1,2,1,1} = zm_spin;	% Z chain
                op.h2term{2,1,1,1} = para.chain{1}.t(1).*sigmaZ./2; op.h2term{2,2,1,1} = zm_spin;
                op.h2term{3,1,1,1} = zm; op.h2term{3,2,1,1} = zm;	% ancilla chain
                op.h2term{4,1,1,1} = zm; op.h2term{4,2,1,1} = zm;
            case para.L
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
% 				[bpz,  ~, ~] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
                if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
                end
                zmx = sparse(size(bpx,1),size(bpx,1));
				zmz = []; %sparse(size(bpz,1),size(bpz,1));
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = zmz;
                op.h2term{1,1,s,1} = zmx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = zmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = zmz; op.h2term{3,2,s,2} = zmz;
                op.h2term{4,1,s,2} = zmz; op.h2term{4,2,s,2} = zmz;
            otherwise
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
% 				[bpz,  ~, ~] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
				if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
				end
				zmz = []; %sparse(size(bpz,1),size(bpz,1));
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = zmz;
                op.h2term{1,1,s,1} = para.chain{1}.t(s).*bpx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = para.chain{1}.t(s).*bmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = zmz; op.h2term{3,2,s,2} = zmz;
                op.h2term{4,1,s,2} = zmz; op.h2term{4,2,s,2} = zmz;							% replace zmz by [] ??
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'SpinBoson3C'
        %%%%%%%%%%%%%%%%%%% Spin-Boson Model - 3-Chain %%%%%%%%%%%%%%%%%%%%%%
		% Not linear, but in multi-chain configuration!
		% Benchmark for multi-chain method! Should be similar to
		%		'SpinBoson2C'
		% Spin is always in chain 1 for backward compatibility
		%
		% working (23/07/15)
		% Created 17/07/15 by F.S.
        switch s
            case 1
				[sigmaX,~,sigmaZ]  = spinop(para.spinbase);
                zm_spin			   = zeros(2);
                op.h1term{1,1}     = - para.hx./2.*sigmaX - para.hz./2.*sigmaZ;
                op.h2term{1,1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{1,2,1,1} = zm_spin;	% X chain
                op.h2term{2,1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{2,2,1,1} = zm_spin;
                op.h2term{3,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{3,2,1,1} = zm_spin;	% Z chain
                op.h2term{4,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{4,2,1,1} = zm_spin;
				op.h2term{5,1,1,1} = para.chain{3}.t(1).*sigmaX./2; op.h2term{5,2,1,1} = zm_spin;	% Z chain
                op.h2term{6,1,1,1} = para.chain{3}.t(1).*sigmaX./2; op.h2term{6,2,1,1} = zm_spin;
            case para.L
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
				[bpz,bmz,nz] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
				[bpy,bmy,ny] = bosonop(para.dk(3,s),para.shift(3,s),para.parity);
                if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
                end
                zmx = sparse(size(bpx,1),size(bpx,1));
				zmz = sparse(size(bpz,1),size(bpz,1));
				zmy = sparse(size(bpy,1),size(bpy,1));
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = para.chain{2}.epsilon(s-1).*nz;
				op.h1term{3,s}     = para.chain{3}.epsilon(s-1).*ny;
                op.h2term{1,1,s,1} = zmx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = zmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = zmz; op.h2term{3,2,s,2} = bmz;
                op.h2term{4,1,s,2} = zmz; op.h2term{4,2,s,2} = bpz;
				op.h2term{5,1,s,3} = zmy; op.h2term{5,2,s,3} = bmy;
                op.h2term{6,1,s,3} = zmy; op.h2term{6,2,s,3} = bpy;
            otherwise
                [bpx,bmx,nx] = bosonop(para.dk(1,s),para.shift(1,s),para.parity);
				[bpz,bmz,nz] = bosonop(para.dk(2,s),para.shift(2,s),para.parity);
				[bpy,bmy,ny] = bosonop(para.dk(3,s),para.shift(3,s),para.parity);
				if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
				end
                op.h1term{1,s}     = para.chain{1}.epsilon(s-1).*nx;
				op.h1term{2,s}     = para.chain{2}.epsilon(s-1).*nz;
				op.h1term{3,s}     = para.chain{3}.epsilon(s-1).*ny;
                op.h2term{1,1,s,1} = para.chain{1}.t(s).*bpx; op.h2term{1,2,s,1} = bmx;
                op.h2term{2,1,s,1} = para.chain{1}.t(s).*bmx; op.h2term{2,2,s,1} = bpx;
                op.h2term{3,1,s,2} = para.chain{2}.t(s).*bpz; op.h2term{3,2,s,2} = bmz;
                op.h2term{4,1,s,2} = para.chain{2}.t(s).*bmz; op.h2term{4,2,s,2} = bpz;
				op.h2term{5,1,s,3} = para.chain{3}.t(s).*bpy; op.h2term{5,2,s,3} = bmy;
                op.h2term{6,1,s,3} = para.chain{3}.t(s).*bmy; op.h2term{6,2,s,3} = bpy;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	case 'SpinBoson5C'
        %%%%%%%%%%%%%%%%%%% Spin-Boson Model - 5-Chain %%%%%%%%%%%%%%%%%%%%%%
		% Not linear, but in multi-chain configuration!
		% Spin is always in chain 1 for backward compatibility
		%
		% working??
		% Created 17/08/15 by F.S.
        switch s
            case 1		% will be pentacene system!
				[sigmaX,~,sigmaZ]   = spinop(para.spinbase);
				zm_spin			    = zeros(2);
				op.h1term{1,1}      = - para.hx./2.*sigmaX - para.hz./2.*sigmaZ;
				op.h2term{1 ,1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{1 ,2,1,1} = zm_spin;	% X chain
				op.h2term{2 ,1,1,1} = para.chain{1}.t(1).*sigmaX./2; op.h2term{2 ,2,1,1} = zm_spin;
				op.h2term{3 ,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{3 ,2,1,1} = zm_spin;	% Z chain
				op.h2term{4 ,1,1,1} = para.chain{2}.t(1).*sigmaZ./2; op.h2term{4 ,2,1,1} = zm_spin;
				op.h2term{5 ,1,1,1} = para.chain{3}.t(1).*sigmaX./2; op.h2term{5 ,2,1,1} = zm_spin;	% Z chain
				op.h2term{6 ,1,1,1} = para.chain{3}.t(1).*sigmaX./2; op.h2term{6 ,2,1,1} = zm_spin;
				op.h2term{7 ,1,1,1} = para.chain{4}.t(1).*sigmaZ./2; op.h2term{7 ,2,1,1} = zm_spin;	% Z chain
				op.h2term{8 ,1,1,1} = para.chain{4}.t(1).*sigmaZ./2; op.h2term{8 ,2,1,1} = zm_spin;
				op.h2term{9 ,1,1,1} = para.chain{5}.t(1).*sigmaX./2; op.h2term{9 ,2,1,1} = zm_spin;	% Z chain
				op.h2term{10,1,1,1} = para.chain{5}.t(1).*sigmaX./2; op.h2term{10,2,1,1} = zm_spin;
            case para.L
				if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
				end
				for i = 1:5			% slow, but easy to modify!
					[bp,bm,n] = bosonop(para.dk(i,s),para.shift(i,s),para.parity);
					zm = sparse(size(bp,1),size(bp,1));
					op.h1term{i,s}		   = para.chain{i}.epsilon(s-1).*n;
					op.h1term{2*i-1,1,s,i} = zm; op.h1term{2*i-1,2,s,i} = bm;
					op.h1term{2*i  ,1,s,i} = zm; op.h1term{2*i  ,2,s,i} = bp;
				end
			otherwise
				if para.parity ~= 'n'
					error('VMPS:genh1h2term_onesite:ParityNotSupported','parity not implemented yet');
				end
				for i = 1:5
					[bp,bm,n] = bosonop(para.dk(i,s),para.shift(i,s),para.parity);
					op.h1term{i,s}		   = para.chain{i}.epsilon(s-1).*n;
					op.h1term{2*i-1,1,s,i} = para.chain{i}.t(s).*bp; op.h1term{2*i-1,2,s,i} = bm;
					op.h1term{2*i  ,1,s,i} = para.chain{i}.t(s).*bm; op.h1term{2*i  ,2,s,i} = bp;
				end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
end

function [H0, H1] = MLSB_Operators(para)
% calculates Energy levels H0 and couplings to bath HI
% define modes and parameters in VMPS_MLSBM
%
% perhaps export to file
switch para.MLSB_mode
    case 1
        %   Define equal spacing Delta between each level:  para.MLSB_Delta
        %   Define couplings System-Bath as vector:         para.MLSB_t
        %   Define Size of ML-System:                       para.MLSBM_Ls
        %   Energy levels will be symmetrically aligned around 0
        %   TODO: next-neighbour couping within system
        assert(length(para.MLSB_Delta) == 1, 'Only one spacing Delta allowed');
        assert(length(para.MLSB_t) == para.MLSB_Ls, 'All couplings t between system and bath must be defined');
        H0 = diag(((para.MLSB_Ls:-1:1)- (para.MLSB_Ls+1)/2).*para.MLSB_Delta);
        H1 = diag(para.MLSB_t);
    case 2
        % Hamiltonian with rotational symmetry. Get from separate function
        % untis in eV
        [H0, H1] = Hamiltonian_PPC(para);
	case 3
		% Hamiltonian of DP-MES
		H0 = zeros(5);
		H1 = zeros(5);
		error('VMPS:genh1h2term_onesite','DP_MES Hamiltonian not yet available.')
    otherwise
end

end
