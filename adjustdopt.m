function [op,para,results,mps,Vmat]=adjustdopt(op,para,results,mps,Vmat)
%Adaptively change $d_{opt}$ so that the smallest sigular value s is below the thredhold svmintol
%Cheng Guo
%17 Feb 2012
% Modified:
%	FS 03/02/2014:	- added abs() for calculation of d_opt_change, as this prevents elimination of positive and negative contributions.


if para.logging													% save d_opt history
		results.d_opt{para.loop-1}=para.d_opt;
end

d_opt_old=para.d_opt;
for s=2:min(para.trustsite(end)+0,para.L) 			% Only modify boson sites
        [Db,d_opt]=size(Vmat{s});
        discarddims=find(results.Vmat_sv{s}<para.svmintol);			% lowest SV below certain threshold --> discard these d_opt dimensions
        if ~isempty(discarddims)
            dispif('Discard Dims',para.logging)
            if para.parity~='n'
                error('Haven not implemented parity code for adaptive d_opt!');
%                 adddim=ceil(maxp*para.dk(s)/2); %This is according to experience,  if maxp=0.5 then this will add 50 new states.
%                 if para.dk(s)+2*adddim>para.dkmax
%                     adddim=(para.dkmax-para.dk(s))/2;
%                 end
%                 para.dk(s)=para.dk(s)+2*adddim;
%                 addmat=zeros(adddim,d_opt);
%                 Vmat{s}=cat(1,addmat,Vmat{s}(1:Db/2,:),addmat,Vmat{s}(Db/2+1:end,:));
%                 Vmat{s}(1:adddim,1:d_opt/2)=1e-15*randn(adddim,d_opt/2); %add a lit perturbation in the new added dims
%                 Vmat{s}(Db/2+adddim+1:Db/2+2*adddim,d_opt/2+1:end)=1e-15*randn(adddim,d_opt/2);
            else if para.d_opt(s)-length(discarddims)>1
                dispif('remove dims in d_opt',para.logging)
                para.d_opt(s)=para.d_opt(s)-length(discarddims);
                Vmat{s}(:,discarddims)=[];
                mps{s}(:,:,discarddims)=[];
                results.Vmat_sv{s}(discarddims)=[];
                else
                    dispif('Done nothing',para.logging)
                end
%                 adddim=ceil(maxp*para.dk(s)); %This is according to experience,  if maxp=0.5 then this will add 50 new states.
%                 if para.dk(s)+adddim>para.dkmax
%                     adddim=para.dkmax-para.dk(s);
%                 end
%                 para.dk(s)=para.dk(s)+adddim;
%                 addmat=zeros(adddim,d_opt);
%                 Vmat{s}=cat(1,addmat,Vmat{s});
            end
        else
            dimincratio=log10(results.Vmat_sv{s}(end)/para.svmaxtol)/2; 	%Imperial relation
            adddim=ceil(para.d_opt(s)*dimincratio); 									% increase dim (only if addim > 0)
			if para.d_opt(s)+adddim >= para.dk(s)
				exceeded = para.d_opt(s)+adddim - para.dk(s) +1;				% By how much dk is too small. +1 as d_opt != dk
				adddim = adddim - exceeded;												% Only increase to: d_opt = dk - 1
				if exceeded > para.increasedk
					para.increasedk = exceeded;											% log highest needed addition
				end
				dispif(['Increase adddim(',num2str(s),') = ',num2str(adddim),'; More needed: ',num2str(exceeded)],para.logging)
			else
				dispif(['Increase adddim(',num2str(s),') = ',num2str(adddim)],para.logging)
			end
            if results.Vmat_sv{s}(end)>para.svmaxtol && para.d_opt(s)+adddim<para.dk(s)  %Expand d_opt
                dispif(['Expanding site ',num2str(s)],para.logging)
                para.d_opt(s)=para.d_opt(s)+adddim;
                addmat=zeros(para.dk(s),adddim);
                Vmat{s}=cat(2,Vmat{s},addmat);
                [a,b,c]=size(mps{s});
                addmat=zeros(a,b,adddim);
                mps{s}=cat(3,mps{s},addmat);
            end
        end
		% expand dk if highest SV < para.expandBelowSV
		if para.useexpand == 1 && max(results.Vmat_sv{s}) < para.expandBelowSV
			maxp = 0.4;							% find a better ratio from SV
			adddim=ceil(maxp*para.dk(s)); %This is according to experience
			 if para.dk(s)+adddim > para.dkmax
				 adddim=para.dkmax-para.dk(s);
			 end
			 para.dk(s) = para.dk(s)+adddim;
			 addmat = zeros(adddim,para.d_opt(s));
			 Vmat{s}=cat(1,addmat,Vmat{s});
			 para.hasexpanded = 1;
			 dispif('Increased dk',para.logging)
		end
end

[op,para]=genh1h2term(para,op);
disp('(para.d_opt-d_opt_old)./para.d_opt=')
(para.d_opt-d_opt_old)./para.d_opt
para.d_opt_change=mean(abs(para.d_opt-d_opt_old)./para.d_opt);

if para.logging													% save d_opt history
		results.d_opt{para.loop}=para.d_opt;
		results.dk{para.loop} = para.dk;
end

end
