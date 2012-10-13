function op=gen_sitej_h1h2(op,para,sitej)

Lambda=para.Lambda;
if sitej>=3 && para.rescaling==1
    op.h1j=Lambda.^(sitej-2).*op.h1term{sitej};
    op.h2j=op.h2term(:,:,sitej);
    for m=1:para.M
    op.h2j{m,1}=op.h2j{m,1}.*Lambda^(sitej-2);
    end
else
    op.h1j=op.h1term{sitej};
    op.h2j=op.h2term(:,:,sitej);
end

end