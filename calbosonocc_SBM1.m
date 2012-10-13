function nx=calbosonocc_SBM1(mps0,Vmat0,para,results)
%Calculate the boson occupation on x chain
%The operater on the spin site is set to zero
n_op=cell(1,para.L);

for j=1:para.L
    if j~=para.spinposition
        [bp,bm,n] = bosonop(para.dk(j),para.shift(j),para.parity);
%         if para.parity=='n'
%             idm=eye(size(n));
%             bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
%         else
%             [bpx,bmx,nx,bpy,bmy,ny]=paritykron(bp,para.bosonparity);
%         end
        n_op{j}=n;
    else
        n_op{para.spinposition}=zeros(para.dk(j));
    end
end
nx =expectation_allsites(n_op,mps0,Vmat0);

end