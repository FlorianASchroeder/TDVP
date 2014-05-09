function [mps,Vmat,para,results] = rightnormA(mps,Vmat,para,results)
% adjust Bond dimension D to optimum
% sweeps rl
%
% Modified:
%	FS 03/02/2014:	- added abs() for calculation of D_change, as this prevents elimination of positive and negative contributions.


if (para.logging == 1) && (para.adjust == 1)
	results.D{para.loop-1} = para.D;
end

D_old=para.D;
N = length(mps);
for i = N:-1:2
    if para.adjust==1
        [mps{i}, U,para,results] = prepare_onesite_truncate(mps{i}, 'rl',para,i,results);		% expand or truncate A in D
    else
        [mps{i}, U] = prepare_onesite(mps{i}, 'rl',para,i);
    end
    mps{i-1} = contracttensors(mps{i-1}, 3, 2, U, 2, 1);
    mps{i-1} = permute(mps{i-1}, [1, 3, 2]);
end
if para.adjust==1
		dispif('(para.D-D_old)./para.D = ',para.logging)
		(para.D-D_old)./para.D
		para.D_change=mean(abs(para.D-D_old)./para.D);
		if para.logging == 1
			results.D{para.loop}=para.D;
		end
end
end
