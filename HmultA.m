function outvec=HmultA(invec, op, Dl, Dr, d, M,parity,nonzeroindex)
% calculates H|MPS> for energy minimization.
% invec & outvec are vectorized forms of the optimized MPS-matrix A
%
% Modified
%	- FS 23/11/14: changed last reshape to numel() since faster
%	- FS 05/08/15: changed from h1j, h2j -> h1jOBB, h2jOBB
if parity~='n'
    inA=zeros(Dl,Dr,d);
    inA(nonzeroindex)=invec;
else
    inA=reshape(invec,[Dl,Dr,d]);
end


out=0;

% collect non-interacting parts of Hamiltonian:
outHl = contracttensors(op.Hleft,2,2,inA,3,1);				% out += Hleft_(l',l) * A_(l,r,n) * 1_(r',r) *1_(n',n)
outHr = contracttensors(inA,3,2,op.Hright,2,2);				% ()_(l',n',r') = A_(l,r,n) * Hright_(r',r) * 1_(l',l) *1_(n',n)
outHr = permute(outHr,[1,3,2]);								% += ()_(l',r',n')
out   = out+outHr+outHl;
if ~isempty(op.h1jOBB)
    out = out+contracttensors(inA,3,3,op.h1jOBB,2,2);
end


% collect interacting parts:
for m=1:M
	outSl = 0; outSr = 0;
	if ~isempty(op.h2jOBB{m,2}) && ~isempty(op.Opleft{m})
		outSl = contracttensors(op.Opleft{m},2,2,inA,3,1);
		outSl = contracttensors(outSl,3,3,op.h2jOBB{m,2},2,2);
	end
	if ~isempty(op.h2jOBB{m,1}) && ~isempty(op.Opright{m})
		outSr = contracttensors(inA,3,2,op.Opright{m},2,2);
		outSr = contracttensors(outSr,3,2,op.h2jOBB{m,1},2,2);
	end
    out = out+outSr+outSl;
end


if parity~='n'
    outvec=out(nonzeroindex);
else
    outvec = reshape(out,[numel(out),1]);
end

end