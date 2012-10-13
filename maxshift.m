function para=maxshift(para)
para.maxshift=para.shift;
for k=1:para.L
    [bp,bm,n]=bosonop(para.dk(k),0,para.parity);
    x=sqrt(2)/2.*(bp+bm);
    x=full(x);
    para.maxshift(k)=max(eig(x));
end
end