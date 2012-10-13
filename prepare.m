function [mps,Vmat,para] = prepare(mps,Vmat,para)

N = length(mps);

for i = 1:N-1
    if para.useVmat==1
        [Vmat{i},V] = prepare_onesiteVmat(Vmat{i},para);
        mps{i} = contracttensors(mps{i},3,3,V,2,2);
    end
    [mps{i}, U] = prepare_onesite(mps{i}, 'lr',para,i);
    mps{i+1} = contracttensors(U,2,2,mps{i+1},3,1);
    para=gennonzeroindex(mps,Vmat,para,i);
    para=gennonzeroindex(mps,Vmat,para,i+1);
end


for i = N:-1:2
%     if para.useVmat==1
%         [Vmat{i},V] = prepare_onesiteVmat(Vmat{i},para);
%         mps{i} = contracttensors(mps{i},3,3,V,2,2);
%     end
    [mps{i}, U] = prepare_onesite(mps{i}, 'rl',para,i);
    mps{i-1} = contracttensors(mps{i-1}, 3, 2, U, 2, 1);
    mps{i-1} = permute(mps{i-1}, [1, 3, 2]);
        para=gennonzeroindex(mps,Vmat,para,i);
        para=gennonzeroindex(mps,Vmat,para,i-1);
end



end