function vNE = vonNeumannEntropy(S)

dim=size(S,1);
s=diag(S);
s=s./sqrt(s'*s);

%maxvNE=log(dim)/log(2.0);

entropy=0;
for x=1:dim
    if(s(x)>1e-16)
        entropy=entropy-s(x).^2.*(log(s(x).^2));
    end
end

vNE=entropy;
