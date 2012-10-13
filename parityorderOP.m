  function newmat=parityorderOP(mat)

%Reordering the boson local operators in to 2*2 block matrix according to
%parity odd and even. The shpae of the block matrix is :
% oo oe
% eo ee
% So operator
% 0 sqrt(3) 0 0
% 0 0 sqrt(2) 0
% 0 0 0 1
% 0 0 0 0
% will become:
% 0 0 sqrt(3) 0
% 0 0 0 1 0
% 0 sqrt(2) 0 0
% 0 0 0 0

[m,n]=size(mat);
assert(m==n);
assert(mod(m,2)==0);
blockdim=m/2;

o=1:2:m;
e=2:2:m;
newmat=mat;
newmat(1:blockdim,1:blockdim)=mat(o,o);
newmat(1:blockdim,blockdim+1:end)=mat(o,e);
newmat(blockdim+1:end,1:blockdim)=mat(e,o);
newmat(blockdim+1:end,blockdim+1:end)=mat(e,e);

end