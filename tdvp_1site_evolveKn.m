function [mps, Vmat] = tdvp_1site_evolveKn(mps,Vmat,para,results,op,sitej,Cn,Hn)
%% Evolves the non-site center C(n) following Haegeman 2014
%   - Only contains C = exp(i K(n) dt/2) C;
%   - Varies for l->r and l<-r; uses para.sweepto
%   - Hn input as (l*r*n,l*r*n) 2D-array
%
% Created by Florian Schroeder @ Cambridge 20/10/2014
switch para.sweepto
    case'r'
        %% Get Dimensions
        BondDimLeft  = size(op.Hlrstorage{sitej},1);
        BondDimRight = results.D{end}(sitej);   % or size(op.Hlrstorage{sitej+1},1), if results.D is not properly written. But should work always since any last change to D is also last entry in that cell array!
        dk           = results.dk{end}(sitej);
        if para.useVmat
            OBBDim      = results.d_opt{end}(sitej);
        else
            OBBDim      = dk;
        end

        %% prepare K(n) from H(n)
        % through contraction with A(t+dt)
        % K(n)_(rl',r',rl,r) = [A*_(l',rl',n') * H(n)_(l',r',n',l,r,n)]_(rl',r',l,r,n) * A_(l,rl,n)
        Hn = reshape(Hn, [BondDimLeft, BondDimRight, OBBDim, BondDimLeft, BondDimRight, OBBDim]);
        Kn = contracttensors(conj(mps{sitej}),3,[1 3], Hn,6,[1 3]);
        Kn = contracttensors(Kn,5,[3 5], mps{sitej},3,[1 3]);                   % K(n)_(rl',r',r,rl)
        Kn = permute(Kn,[1,2,4,3]);                                             % K(n)_(rl',r',rl,r)
        Kn = reshape(Kn,[BondDimRight^2, BondDimRight^2]);                      % K(n)_(rl'r',rlr)

        %% Take and apply Matrix exponential
        % C(n,t) = exp(+ i K(n) dt)_(rl'*r',rl*r) * C(n,t+dt)_(rl*r)
        if size(Kn,1) <= para.tdvp.maxExpMDim
            Cn = expm( 1i .* Kn .* para.tdvp.deltaT./2) *reshape(Cn,[BondDimRight^2,1]);
        else
            Cn = expv(1i*para.tdvp.deltaT./2, Kn,...
                reshape(Cn,[BondDimRight^2,1]),...
                para.tdvp.expvTol, para.tdvp.expvM);
        end
        Cn = reshape(Cn, [BondDimRight, BondDimRight]);
        clear('Kn', 'Hn');

        %% Multiply Ac(n+1) = C(n) * Ar(n+1)
        % set focus on next site A
        % A_(rl,r,n) = C(n)_(rl,l) * A_(l,r,n)
        mps{sitej+1} = contracttensors(Cn,2,2, mps{sitej+1},3,1);
        clear('Cn');
    case 'l'
        %% Get Dimensions
        [OldDimLeft, OldDimRight, OldOBBDim] = size(mps{sitej+1});              % needed to reshape H(n+1)
%         BondDimRight = OldDimLeft;                                              % should always hold! since A|l-r|A

        %% prepare K(n) from H(n+1)
        % through contraction with A(t+dt)
        % H(n) still persistent from previous sweep
        % K(n)_(l',lr',l,lr) = [A*_(lr',r',n') * H(n+1)_(l',r',n',l,r,n)]_(lr',l',l,r,n) * A_(lr,r,n)
        Hn = reshape(Hn, [OldDimLeft, OldDimRight, OldOBBDim, OldDimLeft, OldDimRight, OldOBBDim]);
        Kn = contracttensors(conj(mps{sitej+1}),3,[2 3], Hn,6,[2 3]);       % K(n)_(lr',l',l,r,n)
        Kn = contracttensors(Kn,5,[4 5], mps{sitej+1},3,[2 3]);             % K(n)_(lr',l',l,lr)
        Kn = permute(Kn,[2,1,3,4]);                                         % K(n)_(l',lr',l,lr)
        Kn = reshape(Kn,[OldDimLeft^2, OldDimLeft^2]);                      % K(n)_(l'*lr',l*lr)

        %% Take and apply Matrix exponential
        % C(n,t+dt/2) = exp(+ i K(n) dt/2)_(l'*lr',l*lr) * C(n,t+dt)_(l*lr)
        if size(Kn,1) <= para.tdvp.maxExpMDim
            Cn = expm( 1i .* Kn .* para.tdvp.deltaT./2) *reshape(Cn,[OldDimLeft^2,1]);
        else
            Cn = expv(1i*para.tdvp.deltaT./2, Kn,...
                reshape(Cn,[OldDimLeft^2,1]),...
                para.tdvp.expvTol, para.tdvp.expvM);
        end
        Cn = reshape(Cn, [OldDimLeft, OldDimLeft]);
        clear('Kn', 'Hn');

        %% Multiply Ac(n) = Al(n) * C(n)
        % set focus on next site A
        % A_(l,r,n) = A_(l,r',n) * C(n)_(r',r)
        mps{sitej} = contracttensors(mps{sitej},3,2, Cn,2,1);               % A_(l,n,r)
        mps{sitej} = permute(mps{sitej},[1,3,2]);                           % A_(l,r,n)
        clear('Cn');
end
end