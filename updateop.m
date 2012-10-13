function op=updateop(op,mps,Vmat,sitej,para)
%Update op when finished calculating at site sitej

% if strcmp(hand,'lr')
%     if sitej>=3
%         op.Hleft = op.Hlrstorage{sitej};
%     end
%     op.Hlrstorage{sitej + 1} = updateHleft(op.Hleft, op.h1term{sitej}, op.Opstorage(:,1,sitej), mps{sitej}, Vmat{sitej}, op.h2term(:,2,sitej), mps{sitej}, Vmat{sitej}, M);
%     for m = 1:M
%         op.Opstorage{m,1,sitej+1}=updateCleft([],mps{sitej},Vmat{sitej},op.h2term{m,1,sitej},mps{sitej},Vmat{sitej});
%     end
% else
%     if sitej>=3
%         op.Hright = op.Hlrstorage{sitej + 1};
%     end
%     op.Hlrstorage{sitej} = updateHright(op.Hright, op.h1term{sitej}, op.Opstorage(:,2,sitej+1),mps{sitej},Vmat{sitej}, op.h2term(:,1,sitej), mps{sitej},Vmat{sitej}, M);
%     for m = 1:M
%         h = reshape(hset{m,sitej}, [1, 1, d_opt, d_opt]);
%         op.Opstorage{m, 2, sitej} = updateCright([], mps{sitej},Vmat{sitej}, op.h2term{m,2,sitej}, mps{sitej},Vmat{sitej});
%     end
% end
M=para.M;
if para.sweepto=='r'
    op.Hlrstorage{sitej + 1} = updateHleft(op.Hleft, op.h1j, op.Opleft(:), mps{sitej}, Vmat{sitej}, op.h2j(:,2), mps{sitej}, Vmat{sitej}, M);
    for m = 1:M
        op.Opstorage{m,1,sitej+1}=updateCleft([],mps{sitej},Vmat{sitej},op.h2term{m,1,sitej},mps{sitej},Vmat{sitej});
    end
else
    if sitej>=3 && para.rescaling==1
        op.Hright = op.Hlrstorage{sitej + 1};
    end
    op.Hlrstorage{sitej} = updateHright(op.Hright, op.h1term{sitej}, op.Opstorage(:,2,sitej+1),mps{sitej},Vmat{sitej}, op.h2term(:,1,sitej), mps{sitej},Vmat{sitej}, M);
    for m = 1:M
        op.Opstorage{m, 2, sitej} = updateCright([], mps{sitej},Vmat{sitej}, op.h2term{m,2,sitej}, mps{sitej},Vmat{sitej});
    end
end
end