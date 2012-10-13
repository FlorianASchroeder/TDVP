function [Aj,Vmat,results,para,op]=optimizesite(mps,Vmat,op,para,results,sitej)

optV=1;
while optV
    if para.dk(sitej)>2 &&para.useVmat==1 %Only optizime Vmat{sitej} for the boson sites.
        [Amat,V] = prepare_onesiteAmat(mps{sitej},para,sitej);
        Blaststep = contracttensors(Vmat{sitej}, 2, 2, V, 2, 2);
        [B,E] = minimizeE_onesiteVmat(op, Amat, Blaststep,para);
%         if sitej>=3 && para.rescaling==1
%             results.geoffset(sitej)=(results.geoffset(sitej-1)+results.leftge(sitej))*para.Lambda;
%             E=(E+results.geoffset(sitej))./(para.Lambda.^(sitej-2));
%         end
%         %results.leftge(sitej)
%         fprintf('\nE = %.10g\n', E);
        [Vmat{sitej}, V, results] = prepare_onesiteVmat(B,para,results,sitej);
        Amatlaststep = contracttensors(Amat, 3, 3, V, 2, 2);
    else
        Amatlaststep = mps{sitej}; optV=0;
    end
    [Aj, E] = minimizeE_onesiteA(op, Vmat{sitej}, Amatlaststep,para,sitej);
    if para.loop>1 && para.dk(sitej)>2 && para.useshift==1 && sitej<=para.trustsite(end)
        switch para.model
            case 'SpinDoubleBoson'
                [bp,bm,n] = bosonop(sqrt(para.dk(sitej)),para.shift(sitej),para.parity);
                if para.parity=='n'
                    idm=eye(size(n));
                    bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
                    bpz=kron(idm,bp);bmz=bpz';nz=kron(idm,n);
                else
                    [bpx,bmx,nx,bpz,bmz,nz]=paritykron(bp,para.bosonparity);
                end
                x=sqrt(2)/2*(bpx+bmx);
            otherwise
                [bp,bm,n] = bosonop(para.dk(sitej),para.shift(sitej),para.parity);
                x=sqrt(2)/2*(bp+bm);
        end
        if para.useVmat==1
            temp=contracttensors(x,2,2,Vmat{sitej},2,1);
            temp=contracttensors(conj(Vmat{sitej}),2,1,temp,2,1);
        else
            temp=x;
        end
        temp=contracttensors(Aj,3,3,temp,2,2);
        shift=real(contracttensors(temp,3,[1,2,3],conj(Aj),3,[1,2,3]));
        para.relativeshift(sitej)=abs(shift-para.shift(sitej))/para.maxshift(sitej);
        if  para.relativeshift(sitej)>para.relativeshiftprecision
            para.shift(sitej)=shift;%disp(shift)
            op=update_sitej_h1h2(para,op,sitej);
            mps{sitej}=Aj;
        else
            optV=0;
        end
    else
        optV=0;
    end
end

if sitej>=3 && para.rescaling==1
    results.geoffset(sitej)=(results.geoffset(sitej-1)+results.leftge(sitej))*para.Lambda;
    E=(E+results.geoffset(sitej))./(para.Lambda.^(sitej-2));
end
%results.leftge(sitej)
%fprintf('E = %.10g\n', E);
results.E=E;

end
