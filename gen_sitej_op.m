function op=gen_sitej_op(op,para,sitej,leftge,hand)
%Update op when calculating at site sitej
Lambda=para.Lambda;
op.Opleft = op.Opstorage(:,1,sitej);
op.Opright = op.Opstorage(:,2,sitej+1);
if sitej>=3 && para.rescaling==1
    %op.Hleft =Lambda.^(sitej-2).*(op.Hlrstorage{sitej}-leftge(sitej).*eye(size(op.Hlrstorage{sitej})));
    op.Hleft = Lambda.*(op.Hlrstorage{sitej}-leftge(sitej).*eye(size(op.Hlrstorage{sitej})));
    op.Hright = Lambda.^(sitej-2).*op.Hlrstorage{sitej + 1};
    for m=1:para.M
    op.Opleft{m} = op.Opleft{m}.*Lambda^(sitej-2);
    end
else
    op.Hleft = op.Hlrstorage{sitej};
    op.Hright = op.Hlrstorage{sitej + 1};
end

op=gen_sitej_h1h2(op,para,sitej);

end
