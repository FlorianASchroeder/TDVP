function outvec=HmultA(invec, op, Dl, Dr, d, M,parity,nonzeroindex)

if parity~='n'
    inA=zeros(Dl,Dr,d);
    inA(nonzeroindex)=invec;
else
    inA=reshape(invec,[Dl,Dr,d]);
end


out=0;

outHl=contracttensors(op.Hleft,2,2,inA,3,1);
outHr=contracttensors(inA,3,2,op.Hright,2,2);
outHr=permute(outHr,[1,3,2]);
out=out+outHr+outHl;

if ~isempty(op.h1j)
    out=out+contracttensors(inA,3,3,op.h1j,2,2);
end

for m=1:M
    outSl=contracttensors(op.Opleft{m},2,2,inA,3,1);
    outSl=contracttensors(outSl,3,3,op.h2j{m,2},2,2);

    outSr=contracttensors(inA,3,2,op.Opright{m},2,2);
    outSr=contracttensors(outSr,3,2,op.h2j{m,1},2,2);

    out=out+outSr+outSl;
end


if parity~='n'
    outvec=out(nonzeroindex);
else
    outvec=reshape(out,[Dl*Dr*d,1]);
end

end