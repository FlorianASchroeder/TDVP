function [Aj,Vmat,results,para,op]=optimizesite(mps,Vmat,op,para,results,sitej)
%
% If at spinsite, sitej = 1: Amatlaststep = mps{sitej}; [Aj, E] = minimizeE_onesiteA(op, Vmat{sitej}, Amatlaststep,para,sitej);
% TODO: Improve speed for optV loop by using QR.
%		SVD result needed in lines: 77-79, so perhaps do it before that in a single run only if useshift == 1?
%		    will produce overhead if shift not used as (QR+SVD), but great reduction with loops >=2
%		QR usable in Lines: 13, 22
%
% Commented by Florian Schroeder 13/01/2014
optV=1;
while optV
    if para.dk(sitej)>2 &&para.useVmat==1 % Only use Vmat{j} and optimize for the boson sites.
        [Amat,V] = prepare_onesiteAmat(mps{sitej},para,sitej);							% left-normalize A, as r -> l sweep did right normalize.
        Blaststep = contracttensors(Vmat{sitej}, 2, 2, V, 2, 2);							% set focus on Vmat
        [B,E] = minimizeE_onesiteVmat(op, Amat, Blaststep,para);							% first Energy Optimization
%         if sitej>=3 && para.rescaling==1
%             results.geoffset(sitej)=(results.geoffset(sitej-1)+results.leftge(sitej))*para.Lambda;
%             E=(E+results.geoffset(sitej))./(para.Lambda.^(sitej-2));
%         end
%         %results.leftge(sitej)
%         fprintf('\nE = %.10g\n', E);
        [Vmat{sitej}, V, results] = prepare_onesiteVmat(B,para,results,sitej);
        Amatlaststep = contracttensors(Amat, 3, 3, V, 2, 2);
    else
        Amatlaststep = mps{sitej}; optV=0;
    end
    [Aj, E] = minimizeE_onesiteA(op, Vmat{sitej}, Amatlaststep,para,sitej);

% Shifting of bosons if trustsite > 0;  modify condition for sitej = 1, dk(1) = 4; perhaps: para.dk(sitej) > para.dk(para.spinposition)
% Necessary: dk(j) always > spin dimension
    if para.loop>1 && para.dk(sitej)>2 && para.useshift==1
		if para.useFloShift == 1 && ~(para.trustsite(end)>0)						% stop loops if criterium not fulfilled
			optV = 0;
%		elseif para.useFloShift2 == 1 && ~(min(cellfun(@(x) x(1,1), results.Vmat_sv(2:end))) < para.FloShift2minSV)	% cellfun takes all maximum SV of each site into an array
		elseif para.useFloShift2 == 1 && ~(min(results.Vmat_svLog{end}) < para.FloShift2minSV)	% results.Vmat_svLog{end} contains max SV from last sweep
			optV = 0;
%		elseif para.useFloShift3 == 1 && (~(mod(para.loop,para.FloShift3loops)==0) || ~(results.Vmat_svLog{end}(sitej-1) < (1.+para.FloShift3minMaxSV)/2.))		%FloShift3minMaxSV from last sweep. (sitej-1) as spinsite not included
		elseif para.useFloShift3 == 1 && (~(mod(para.loop,para.FloShift3loops)==0) || ~((results.Vmat_svLog{end}(sitej-1) < (2.+para.FloShift3minMaxSV)/3.) || (sitej-1 < find(results.Vmat_svLog{end} <= para.FloShift3minMaxSV,1) && results.Vmat_svLog{end}(sitej-1) ~= 1)))
		% Apply shift if: every nth loop; (site has had a smaller SV than criterium) or (site is below worst site and SV ~= 1)
			optV = 0;
		elseif para.useChengShift == 1 && ~(sitej<=para.trustsite(end))		% -- original condition from Cheng.
			optV = 0;
		elseif ~para.useEveryShift && ~para.useChengShift && ~para.useFloShift && ~para.useFloShift2 && ~para.useFloShift3	% wrong definition
			optV = 0;
			disp('Define a shift method!');
	% && results.Eerror < sqrt(para.precision)  ---  inserted energy criterium myself!
		else																							% apply shift
			switch para.model
				case 'SpinDoubleBoson'
					[bp,bm,n] = bosonop(sqrt(para.dk(sitej)),para.shift(sitej),para.parity);
					if para.parity=='n'
						idm=eye(size(n));
						bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
						bpz=kron(idm,bp);bmz=bpz';nz=kron(idm,n);
					else
						[bpx,bmx,nx,bpz,bmz,nz]=paritykron(bp,para.bosonparity);
					end
					x=sqrt(2)/2*(bpx+bmx);
				case '2SpinPhononModel'
					[bp,bm,n] = bosonop(sqrt(para.dk(sitej)),para.shift(sitej),para.parity);
					if para.parity=='n'
						idm=eye(size(n));
						bpr=kron(bp,idm);   bmr=bpr';   nr=kron(n,idm);		% right chain
						bpl=kron(idm,bp);   bml=bpl';   nl=kron(idm,n);		% left chain
					else
						[bpr,bmr,nr,bpl,bml,nl]=paritykron(bp,para.bosonparity);
					end
					x = sqrt(2)/2*(bpr+bmr);								% why only evaluate it for right part?
				otherwise
					[bp,bm,n] = bosonop(para.dk(sitej),para.shift(sitej),para.parity);
					x=sqrt(2)/2*(bp+bm);
			end
			if para.useVmat==1
				temp=contracttensors(x,2,2,Vmat{sitej},2,1);
				temp=contracttensors(conj(Vmat{sitej}),2,1,temp,2,1);
			else
				temp=x;
			end
			temp=contracttensors(Aj,3,3,temp,2,2);
			shift=real(contracttensors(temp,3,[1,2,3],conj(Aj),3,[1,2,3]));					% calculate shift after current optimization
			para.relativeshift(sitej)=abs(shift-para.shift(sitej))/para.maxshift(sitej);
			if  para.relativeshift(sitej)>para.relativeshiftprecision									% if relevant shift, update and rerun optimization
				para.shift(sitej)=shift;
				para.shifted = 1;
				op=update_sitej_h1h2(para,op,sitej);		% shift boson operators
				mps{sitej}=Aj;								% store Aj in mps to use in next loop
			else
				optV=0;
			end
		end
    else
        optV=0;
    end
end

if sitej>=3 && para.rescaling==1
    results.geoffset(sitej)=(results.geoffset(sitej-1)+results.leftge(sitej))*para.Lambda;
    E=(E+results.geoffset(sitej))./(para.Lambda.^(sitej-2));
end
%results.leftge(sitej)
%fprintf('E = %.10g\n', E);
results.E=E;

end
