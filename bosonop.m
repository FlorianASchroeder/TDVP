function [bp,bm,n]=bosonop(dim,shift,parity)
%Generate local boson operators
bp=sparse(dim,dim);
n=1;
for i=dim-1:-1:1
    bp(i,i+1)=sqrt(n);
    n=n+1;
end

if parity~='n'
    bp=parityorderOP(bp);
end

bm=bp';
n=bp*bm;

if shift~=0
    bp=full(bp);bm=full(bm);n=full(n);
    x=sqrt(2)/2.*(bp+bm);
    iden=eye(dim);
    bp=bp+sqrt(2)/2.*shift*iden;
    bm=bm+sqrt(2)/2.*shift*iden;
    n=n+shift.*x+shift^2/2.*iden;
end

end
