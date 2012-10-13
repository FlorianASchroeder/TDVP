function bosonshift=calbosonshift_SBM1(mps0,Vmat0,para,results)
%Calculate boson shift x, x^2, var(x)
%The operater on the spin site is set to zero

x_opx=cell(1,para.L);
x2_opx=cell(1,para.L);
for j=1:para.L
    if j~=para.spinposition
        [bp,bm,n] = bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
        idm=eye(size(n));
        bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
        x_opx{j}=sqrt(2)/2.*(bpx+bmx);
        x2_opx{j}=x_opx{j}*x_opx{j};
    else
        x_opx{j}=zeros(para.dk(j));
        x2_opx{j}=zeros(para.dk(j));
    end
end
bosonshift.x =expectation_allsites(x_opx,mps0,Vmat0);
bosonshift.xsqure = expectation_allsites(x2_opx,mps0,Vmat0);
bosonshift.xvariant = sqrt(bosonshift.xsqure-bosonshift.x.^2);
bosonshift.xerror=mean(abs(para.shift-bosonshift.x));
end