function [mps,Vmat,para,results] = rightnormA(mps,Vmat,para,results)

D_old=para.D;
N = length(mps);
for i = N:-1:2
    if para.adjust==1
        [mps{i}, U,para,results] = prepare_onesite_truncate(mps{i}, 'rl',para,i,results);
    else
        [mps{i}, U] = prepare_onesite(mps{i}, 'rl',para,i);
    end
    mps{i-1} = contracttensors(mps{i-1}, 3, 2, U, 2, 1);
    mps{i-1} = permute(mps{i-1}, [1, 3, 2]);
end
if para.adjust==1
(para.D-D_old)./para.D
para.D_change=mean((para.D-D_old)./para.D);
end
end
