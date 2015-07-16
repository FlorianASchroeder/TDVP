function op=gen_sitej_h1h2(op,para,sitej)
% Take op.h1term{j} -> op.h1j and op.h2term{:,:,j} -> op.h2j
% Apply rescale if wanted! (para.rescale == 1)
%
% Commented by Florian Schroeder 13/01/2014
%                                20/10/2014
% Modified
%	FS 16/07/2015: - added Multi-chain support at OBB level (perhaps not
%						used)

if para.nChains == 1
	op.h1j = op.h1term{sitej};
else
	op.h1j = op.h1term(:,sitej);
end
op.h2j = op.h2term(:,:,sitej,:);

% apply Rescaling to bosonic chain
if sitej>=3 && para.rescaling==1
	assert(para.nChains == 1, 'rescaling not implemented in gen_sitej_h1h2.m yet');
    Lambda=para.Lambda;
    op.h1j=Lambda.^(sitej-2).*op.h1j;		% old: op.h1j=Lambda.^(sitej-2).*op.h1term{sitej};

    for m=1:para.M
        % apply rescaling only to h2term(m,1)?? perhaps since otherwise
        % mixed scaling in interaction terms?
    	op.h2j{m,1}=op.h2j{m,1}.*Lambda^(sitej-2);
    end
end

end