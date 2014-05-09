function op=genh1h2term_onesite(para,op,s)
%Define the hamiltonian
% op.h2term{i,j,k}:
%   k is site number
%   i is number of term in sum to address:  t1*a*b + t1*b*a
%                       sum:                    1       2
%   j is position nr in each hopping term:  t*a*b
%       t1*a*b:     para.t(1)*op.h2term{1,1,1}*op.h2term{1,2,2}
%       t1*b*a:     para.t(1)*op.h2term{2,1,1}*op.h2term{2,2,2}
%example: op.h2term{1,1,1} is with op.h2term{1,2,2}
% interaction terms don't seem to be normal ordered...

% If Model changed, modify also:
%   VMPS_SBM1.m:    para.dk(1) gives dimension of first site;
%   calspin.m:      sx,sy,sz defines spin measure operator. Change if dim(site1) changes!
switch para.model
    case 'SpinBoson'
        %%%%%%%%%%%%%%%%%%%Spin-boson Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch s
            case 1                                                  % first chain pos = all spin sites!
                [sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);       % gives XYZ operators with respect to specific main base
                zm_spin=zeros(2);
                op.h1term{1}=-para.hx./2.*sigmaX-para.hz./2.*sigmaZ;
                op.h2term{1,1,1} = para.t(1).*sigmaZ; op.h2term{1,2,1} = zm_spin;		% was sigmaX
                op.h2term{2,1,1} = para.t(1).*sigmaZ; op.h2term{2,2,1} = zm_spin;		% was sigmaX
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'SpinDoubleBoson'
        %%%%%%%%%%%%%%%%%%%Spin doulbe boson Model One Chain (to each side)%%%%%%%%%%%%%%%%%%%%%%
        switch s
            case 1
                [sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
                zm_spin=zeros(2);
                %assert(para.Delta==0);
                op.h1term{1}=-para.hx./2.*sigmaX-para.hy./2.*sigmaY-para.hz./2.*sigmaZ;
                op.h2term{1,1,1} = para.t(1).*sigmaX; op.h2term{1,2,1} = zm_spin;
                op.h2term{2,1,1} = para.t(1).*sigmaX; op.h2term{2,2,1} = zm_spin;   %to the right
                op.h2term{3,1,1} = para.t(1).*sigmaY; op.h2term{3,2,1} = zm_spin;
                op.h2term{4,1,1} = para.t(1).*sigmaY; op.h2term{4,2,1} = zm_spin;   %to the left
            case para.L
                [bp,bm,n] = bosonop(sqrt(para.dk(para.L)),para.shift(para.L),para.parity);
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
                [bp,bm,n] = bosonop(sqrt(para.dk(s)),para.shift(s),para.parity);
                if para.parity=='n'
                    idm=eye(size(n));
                    bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
                    bpy=kron(idm,bp);bmy=bpy';ny=kron(idm,n);
                else
                    [bpx,bmx,nx,bpy,bmy,ny]=paritykron(bp,para.bosonparity);
                end
                zm=sparse(size(bpx,1),size(bpx,1));
                op.h1term{s}=para.epsilon(s-1).*nx+para.epsilon(s-1).*ny;
                op.h2term{1,1,s} = para.t(s).*bpx; op.h2term{1,2,s} = bmx;
                op.h2term{2,1,s} = para.t(s).*bmx; op.h2term{2,2,s} = bpx;
                op.h2term{3,1,s} = para.t(s).*bpy; op.h2term{3,2,s} = bmy;
                op.h2term{4,1,s} = para.t(s).*bmy; op.h2term{4,2,s} = bpy;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'SpinDoulbeBosonFolded' %%Two chain case. Doesn't work any more.
        %%%%%%%%%%%%%%%%%%Spin Doulbe Boson Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                op.h2term{1,1,s} = para.tz(1).*sigmaZ; op.h2term{1,2,s} = para.ty(1).*sigmaY;
                op.h2term{2,1,s} = para.tz(1).*sigmaZ; op.h2term{2,2,s} = para.ty(1).*sigmaY;
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
                op.h2term{1,1,1} = para.t(1).*(bp+bm); op.h2term{1,2,1} = zm;
                op.h2term{2,1,1} = para.t(1).*(bp+bm); op.h2term{2,2,1} = zm;
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
        %%%%%%%%%%%%%%%%%%%2-State Phonon Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % copied from SpinDoubleBoson. This is indeed the folded model, as the bosonic sites have kron()
        % somehow this looks very much like the folded model because of the kron() product in the boson operators
        switch s
            case 1      % if I want to include two exciton sites with chain-ex1-ex2-chain I need to have ex1,ex2 in site 1 combined!
                        % start off with only one site and extend it for later!
                        % possibly wrong: COUPLING in Prior, Chin 10: parallel to energy splitting! so sigmaZ
                [sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
                sigmaPlus = (sigmaX+i.*sigmaY)./2;
                sigmaMinus = (sigmaX - i.*sigmaY)./2;
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
                [bp,bm,n] = bosonop(sqrt(para.dk(para.L)),para.shift(para.L),para.parity);
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
                [bp,bm,n] = bosonop(sqrt(para.dk(s)),para.shift(s),para.parity);
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

    case '2SpinPhononModelfolded'
        %%%%%%%%%%%%%%%%%%2 Spin Phonon model folded%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2 chains folded into one single chain of dimension^2 for symmetrical coupling
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
                op.h2term{1,1,s} = para.t(1).*sigmaZ; op.h2term{1,2,s} = para.t(1).*sigmaZ;     %use same t in both chains
                op.h2term{2,1,s} = para.t(1).*sigmaZ; op.h2term{2,2,s} = para.t(1).*sigmaZ;
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
    case 'MLSpinBoson'
        %%%%%%%%%%%%%%%%%%%Multi-Level Spin-boson Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for use of this model:
        %   para.hx: de
        switch s
            case 1                                                  % first chain pos = all spin sites!
                [HS0, HSI] = MLSB_Operators(para);                  % uses local function
                zm_spin=zeros(para.MLSB_Ls);
                op.h1term{1}= HS0;
                op.h2term{1,1,1} = para.t(1).*HSI; op.h2term{1,2,1} = zm_spin;		% was sigmaX
                op.h2term{2,1,1} = para.t(1).*HSI; op.h2term{2,2,1} = zm_spin;		% was sigmaX
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

end
end

function [H0, HI] = MLSB_Operators(para)
% calculates Energy levels H0 and couplings to bath HI
% define modes and parameters in VMPS_MLSBM
%
% perhaps export to file
switch para.MLSB_mode
    case 1
        assert(length(para.MLSB_Delta) == 1, 'Only one spacing Delta allowed');
        assert(length(para.MLSB_t) == para.MLSB_Ls, 'All couplings t between system and bath must be defined');
        H0 = diag(((para.MLSB_Ls:-1:1)- (para.MLSB_Ls+1)/2).*para.MLSB_Delta);
        HI = diag(para.MLSB_t);
    otherwise
end

end
