function op=gen_sitej_h1h2(op,para,sitej)
% Take op.h1term{j} -> op.h1j and op.h2term{:,:,j} -> op.h2j
% Apply rescale if wanted! (para.rescale == 1)
%
% Commented by Florian Schroeder 13/01/2014

op.h1j = op.h1term{sitej};
op.h2j = op.h2term(:,:,sitej);

% apply Rescaling to bosonic chain
if sitej>=3 && para.rescaling==1
    Lambda=para.Lambda;
    op.h1j=Lambda.^(sitej-2).*op.h1j;		% old: op.h1j=Lambda.^(sitej-2).*op.h1term{sitej};

    for m=1:para.M
	op.h2j{m,1}=op.h2j{m,1}.*Lambda^(sitej-2);
    end
end

end