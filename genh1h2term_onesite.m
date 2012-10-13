function op=genh1h2term_onesite(para,op,s)
%Define the hamiltonian
switch para.model
    case 'SpinBoson'
        %%%%%%%%%%%%%%%%%%%Spin-boson Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch s
            case 1
                [sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
                zm_spin=zeros(2);
                op.h1term{1}=-para.hx./2.*sigmaX-para.hz./2.*sigmaZ;
                op.h2term{1,1,1} = para.t(1).*sigmaX; op.h2term{1,2,1} = zm_spin;
                op.h2term{2,1,1} = para.t(1).*sigmaX; op.h2term{2,2,1} = zm_spin;
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

    case 'SpinDoubleBoson'
        %%%%%%%%%%%%%%%%%%%Spin doulbe boson Model One Chain%%%%%%%%%%%%%%%%%%%%%%
        switch s
            case 1
                [sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
                zm_spin=zeros(2);
                %assert(para.Delta==0);
                op.h1term{1}=-para.hx./2.*sigmaX-para.hy./2.*sigmaY-para.hz./2.*sigmaZ;
                op.h2term{1,1,1} = para.t(1).*sigmaX; op.h2term{1,2,1} = zm_spin;
                op.h2term{2,1,1} = para.t(1).*sigmaX; op.h2term{2,2,1} = zm_spin;
                op.h2term{3,1,1} = para.t(1).*sigmaY; op.h2term{3,2,1} = zm_spin;
                op.h2term{4,1,1} = para.t(1).*sigmaY; op.h2term{4,2,1} = zm_spin;
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
            case para.spinposition
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
end
end
