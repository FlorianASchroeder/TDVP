function outvec=HmultVmat(invec,op,dkj,d_opt,M,parity)
% invec = inB = optimising Vmat OBB matrix.
%
% Modified:
%	- 20/12/14 FS: faster matrix product for 2-2 contractions
if parity~='n'
    inB=zeros(dkj,d_opt);
    inB(1:dkj/2,1:d_opt/2)=reshape(invec(1:dkj*d_opt/4),dkj/2,d_opt/2);
    inB(dkj/2+1:end,d_opt/2+1:end)=reshape(invec(dkj*d_opt/4+1:end),dkj/2,d_opt/2);
else
    inB=reshape(invec,[dkj,d_opt]);
end

out=0;

outHl = inB * op.HlOPB.';
% outHl = contracttensors(inB,2,2,op.HlOPB,2,2);
outHr = inB * op.HrOPB.';
% outHr = contracttensors(inB,2,2,op.HrOPB,2,2);

out=out+outHl+outHr;

if ~isempty(op.h1j)
	out = out + (op.h1j * inB);
% 	out = out + contracttensors(op.h1j,2,2,inB,2,1);
end

for m=1:M
% 	outSl=contracttensors(op.OpleftOPB{m},2,2,inB,2,2);
% 	outSl=contracttensors(op.h2j{m,2},2,2,outSl,2,2);
	outSl = op.h2j{m,2} * (inB * op.OpleftOPB{m}.');

% 	outSr=contracttensors(inB,2,2,op.OprightOPB{m},2,2);
% 	outSr=contracttensors(op.h2j{m,1},2,2,outSr,2,1);
    outSr = op.h2j{m,1} * (inB * op.OprightOPB{m}.');

    out=out+outSr+outSl;
end

if parity~='n'
    outvec=vertcat(reshape(out(1:dkj/2,1:d_opt/2),dkj*d_opt/4,1),reshape(out(dkj/2+1:end,d_opt/2+1:end),dkj*d_opt/4,1));
else
    outvec=reshape(out,[numel(out),1]);
end

end