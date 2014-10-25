function [op,para,results,mps,Vmat] = adjustdopt(op,para,results,mps,Vmat)
% Adaptively change $d_{opt}$ so that the smallest sigular value s is below the thredhold svmintol
%   Either truncate or expand d_opt per sweep!
% Also change dk using schemes:
%   - para.useDkExpand1
%   - para.useDkExpand2
%
%Cheng Guo
%17 Feb 2012
% Modified:
%	FS 03/02/2014:	- added abs() for calculation of d_opt_change, as this prevents elimination of positive and negative contributions.
%   FS 30/05/2014:  - added that always one sv < svmaxtol is kept, to prevent oscillating expansion / deletion
%                   - added para.d_opt_min as parameter.
%                   - added if-statements to check dimincratio
%                   - Expand dk if tail of Vmat indicates, that expansion is needed.
%   FS 25/10/2014:  - Dk expansion for para.useVmat=0 corrected (eye for
%                     Vmat)
%
%   TODO: replace for by while loop and remove limitation by trustsite, as this adds artificial step.


if para.logging													% save d_opt history
		results.d_opt{para.loop-1} = para.d_opt;
end

d_opt_old = para.d_opt;

for s = 2:min(para.trustsite(end)+0,para.L) 			% Only modify boson sites
    %% Adjust d_opt OBB dimensions
    [Db,d_opt]      = size(Vmat{s});
    if para.useVmat
        para.increasedk = 0;
        discarddims     = find(results.Vmat_sv{s}<para.svmintol);       % lowest SV below certain threshold --> discard these d_opt dimensions

        if ~isempty(discarddims)
            %% Truncate d_opt
            dispif('Discard Dims in d_opt',para.logging)
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
    %            else if para.d_opt(s)-length(discarddims) >= para.d_opt_min            %old 30/05/2014
            else
                % parity == 'n'
                if results.Vmat_sv{s}(discarddims(1)-1) > para.svmaxtol         % if next highest not-discarded element too large (causing expansion in next sweep)
                    discarddims = discarddims(2:end);                           % always keep one element < svmaxtol
                end
                % if Dif > 0: would remove too many dims; if Dif < 0, no
                % problem so set Dif = 0
                difference = para.d_opt_min + length(discarddims) - para.d_opt(s);
                if difference <= 0      % discarddims does not violate d_opt_min
                    difference = 0;
                end                     % else: discarddims would remove too much by amount = diff;
                dispif('remove dims in d_opt',para.logging)
                para.d_opt(s)=para.d_opt(s)-length(discarddims-difference);
                Vmat{s}(:,discarddims(difference+1:end))=[];                          % only cut the end
                mps{s}(:,:,discarddims(difference+1:end))=[];
                results.Vmat_sv{s}(discarddims(difference+1:end))=[];
    %                else                                                       % Don't need this case anymore.
    %                    dispif('Done nothing',para.logging)
    %                end
    %                 adddim=ceil(maxp*para.dk(s)); %This is according to experience,  if maxp=0.5 then this will add 50 new states.
    %                 if para.dk(s)+adddim>para.dkmax
    %                     adddim=para.dkmax-para.dk(s);
    %                 end
    %                 para.dk(s)=para.dk(s)+adddim;
    %                 addmat=zeros(adddim,d_opt);
    %                 Vmat{s}=cat(1,addmat,Vmat{s});
            end
        else
            %% Expand d_opt.
            % from here: no SV < para.svmintol
            dimincratio=log10(results.Vmat_sv{s}(end)/para.svmaxtol)/2; 	%Imperial relation
            if dimincratio > 0
                adddim=ceil(para.d_opt(s)*dimincratio); 									% increase dim (only if addim > 0)
                if para.d_opt(s)+adddim >= para.dk(s)
                    exceeded = para.d_opt(s)+adddim - para.dk(s) +1;				% By how much dk is too small. +1 as d_opt != dk
                    adddim = adddim - exceeded;												% Only increase to: d_opt = dk - 1
                    if exceeded > para.increasedk
                        para.increasedk = exceeded;											% log highest needed addition
                    end
                    dispif(['Increase dopt adddim(',num2str(s),') = ',num2str(adddim),'; More needed: ',num2str(exceeded)],para.logging)
                else
                    dispif(['Increase dopt adddim(',num2str(s),') = ',num2str(adddim)],para.logging)
                end
                if results.Vmat_sv{s}(end)>para.svmaxtol && para.d_opt(s)+adddim<para.dk(s)  % Expand d_opt; this if could be removed.
                    dispif(['Expanding site ',num2str(s)],para.logging)
                    para.d_opt(s)=para.d_opt(s)+adddim;
                    addmat=zeros(para.dk(s),adddim);
                    Vmat{s}=cat(2,Vmat{s},addmat);
                    [a,b,c]=size(mps{s});
                    addmat=zeros(a,b,adddim);
                    mps{s}=cat(3,mps{s},addmat);
                end
            %else: svmaxtol > smallest SV > svmintol --> good, nothing to do
            end
        end
    end
    %% Expand dk
    if para.useexpand == 1
        %% Estimate amount of expansion
        if para.useDkExpand1 == 1 && para.useVmat && max(results.Vmat_sv{s}) < para.expandBelowSV
            % Expand dk if highest SV < para.expandBelowSV
            dimincratio = 0.4;							% find a better ratio from SV
        elseif para.useDkExpand2 == 1 && para.useVmat
            % Expand dk according to Vmat SV-tail analysis
            % vector of length dk:
            % calculate the sum over the weighted Vmat contributions. This tells which states are occupied / needed
            % take min, because dimensions might have changed above.
            occupation = sum(abs(real(Vmat{1,s}(:,1:min(d_opt,size(Vmat{s},2)))*diag(results.Vmat_sv{1,s}(1:min(d_opt,size(Vmat{s},2)))))),2);

            occTail = occupation(1:ceil(para.dkEx2_tail*length(occupation)));                                   % high energy tail
            if std(log10(occTail)) > para.dkEx2_maxDev && abs(mean(log10(occTail))) < para.svmintol^2           % measure order of fluctuations on tail. should be small for flat tail.
                dimincratio = 1/mean(abs(diff(log10(occTail))));                                                % the steeper the less increase needed.
                %dimincratio = log10(std(occTail)/para.dkEx2_maxDev)/3;                                         % similar to imperial relation, slightly modified (estimated). always positive.
            elseif mean(occTail) > para.svmintol                                                                % if absolute value is too high
                dimincratio = 1/mean(abs(diff(log10(occTail))));                                                % some increase
            else
                dimincratio = 0;                        % no increase
            end
        else
            dimincratio = 0;
        end
        dimincratio = max(dimincratio,para.increasedk/para.dk(s));      % if d_opt needs more dk, this assures it.

        %% Apply the expansion
        if dimincratio > 0
            adddim=ceil(dimincratio*para.dk(s));                        % from here: copied from above
            if para.dk(s)+adddim > para.dkmax
                adddim=para.dkmax-para.dk(s);
            end
            para.dk(s) = para.dk(s)+adddim;             % operators will be expanded in genh1h2term?
            if para.useVmat
                addmat = zeros(adddim,para.d_opt(s));
                Vmat{s}=cat(1,addmat,Vmat{s});
            else
                Vmat{s} = eye(para.dk(s));              % TODO: might need changing to eye(para.dk(s), para.d_opt(s)); (see createrandomVmat.m)
            end
            para.hasexpanded = 1;
            dispif('Increased dk',para.logging)
        end
    end
end

% Expand operators
[op,para] = genh1h2term(para,op);

disp('(para.d_opt-d_opt_old)./para.d_opt=')
(para.d_opt-d_opt_old)./para.d_opt
para.d_opt_change=mean(abs(para.d_opt-d_opt_old)./para.d_opt);

if para.logging													% save d_opt history
		results.d_opt{para.loop}=para.d_opt;
		results.dk{para.loop} = para.dk;
end

end
