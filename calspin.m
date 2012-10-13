function spin=calspin(mps,Vmat,para,results)
%Calculate the spin expectation value

N=para.L;
assert(N==length(mps) && N==length(Vmat));

ndset=cell(1,N);

for j=1:N
    ndset{1,j}=eye(size(Vmat{j},1));
end
sx=ndset;sy=sx;sz=sy;

[sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
sx{para.spinposition}=sigmaX;
sy{para.spinposition}=sigmaY;
sz{para.spinposition}=sigmaZ;
spin.sx=expectationvalue(sx,mps,Vmat,mps,Vmat);
spin.sy=expectationvalue(sy,mps,Vmat,mps,Vmat);
spin.sz=expectationvalue(sz,mps,Vmat,mps,Vmat);
spin.sx=real(spin.sx);
spin.sy=real(spin.sy);
spin.sz=real(spin.sz);

end
