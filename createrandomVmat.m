function [Vmat] = createrandomVmat(para)
% Creates random V-matrix used for optimal boson basis.
% Dimensions per site: dk x d_opt
% Modified
%	FS 16/07/2015: - use prod(para.dk(:,I)) to allow multiple chains at OBB level

Vmat = cell(1, para.L);

%Vmat{1}=[1 0;0 1];

for i = 1:para.L
    if para.useVmat==1 && prod(i~= para.spinposition)				% ready for array in spinposition
        if para.parity~='n'
			assert(para.nChains == 1, 'only single chains allowed with this parity');
            halfdkdim=para.dk(i)/2;
            halfddim=para.d_opt(i)/2;
            Vmat{i} = zeros(para.dk(i), para.d_opt(i));
            Vmat{i}(1:halfdkdim,1:halfddim)=randn(halfdkdim,halfddim)./sqrt(halfdkdim*halfddim);
            Vmat{i}(halfdkdim+1:end,halfddim+1:end)=randn(halfdkdim,halfddim)./sqrt(halfdkdim*halfddim);
        else
            Vmat{i} = randn(prod(para.dk(:,i)), para.d_opt(i))./sqrt(prod(para.dk(:,i))*para.d_opt(i));
        end
    else
        Vmat{i} = eye(para.dk(i), para.d_opt(i));                   % why still para.d_opt??
    end
end
end
