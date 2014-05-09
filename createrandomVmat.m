function [Vmat] = createrandomVmat(para)
% Creates random V-matrix used for optimal boson basis.
% Dimensions per site: dk x d_opt
Vmat = cell(1, para.L);

%Vmat{1}=[1 0;0 1];

for i = 1:para.L
    if para.useVmat==1 && i~=1
        if para.parity~='n'
            halfdkdim=para.dk(i)/2;
            halfddim=para.d_opt(i)/2;
            Vmat{i} = zeros(para.dk(i), para.d_opt(i));
            Vmat{i}(1:halfdkdim,1:halfddim)=randn(halfdkdim,halfddim)./sqrt(halfdkdim*halfddim);
            Vmat{i}(halfdkdim+1:end,halfddim+1:end)=randn(halfdkdim,halfddim)./sqrt(halfdkdim*halfddim);
        else
            Vmat{i} = randn(para.dk(i), para.d_opt(i))./sqrt(para.dk(i)*para.d_opt(i));
        end
    else
        Vmat{i} = eye(para.dk(i), para.d_opt(i));
    end
end
end
