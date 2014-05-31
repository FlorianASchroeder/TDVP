function [mps,Vmat,para,results,op] = minimizeE(op,para)
% Ground state calculation
% Modified:
%	FS 30/05/2014: - First verion of a energy flowdiagram added. See Cheng 2012 for more information.
%
% Commented ba Florian Schroeder 03/02/2014

randn('state', 0)
L=para.L;
M = size(op.h2term,1);

if para.resume==1 && para.savedexist==1
    [Vmat,mps,loop,para,results,op]=loadsaved(para);
else
    loop=1;
    Vmat = createrandomVmat(para);
    mps = createrandommps(para);
    %Preassign space for the results structure.
    results=initresults(para);
end

para=gennonzeroindex(mps,Vmat,para);
para

[mps,Vmat,para] = prepare(mps,Vmat,para);
% storage-initialization
[op] = initstorage(mps, Vmat, op,para);

%The relative change of d_opt and D, 1 means 100%.
para.d_opt_change=1;
para.D_change=1;

while loop<=para.loopmax;
    fprintf('\nloop = %d\n',loop);
    para.loop=loop;
	para.shifted = 0;													% to log if bosons were shifted
    results.Evalues = [];
    results.flowdiag{loop} = [];
    % ********* cycle 1: j ? j + 1 (from 1 to L - 1)*********
    for j = 1:L
        fprintf('%g-', j); para.sweepto='r';
        op=gen_sitej_op(op,para,j,results.leftge,'lr');					% take Site h1 & h2 Operators apply rescaling to Hleft, Hright, Opleft ...???
        [Aj,Vmat,results,para,op]=optimizesite(mps,Vmat,op,para,results,j);
        if j~=L
            [mps{j}, U, para,results] = prepare_onesite(Aj, 'lr',para,j,results);
            mps{j+1} = contracttensors(U,2,2,mps{j+1},3,1);
        else
            mps{j}=Aj;
        end
        op=updateop(op,mps,Vmat,j,para);						% calls updateHleft, updateCleft to update for next sweep
		Hleft=op.Hlrstorage{j+1};
		eigvalues=sort(eig((Hleft'+Hleft)/2));
        % log the flowdiagram
        if para.logging && j~=L                                                         % Last site has different length in EV
            if j~=1            % check for different dimensions to apply vertcat
                dimCat = min(size(results.flowdiag{loop},2),size(eigvalues,1));
                results.flowdiag{loop} = [results.flowdiag{loop}(:,1:dimCat); eigvalues(1:dimCat)'];
            else
                results.flowdiag{loop} = [results.flowdiag{loop}; eigvalues'];
            end
        end
		results.leftge(j+1)=eigvalues(1);
        results.Evalues = [results.Evalues, results.E];		% results.E from optimizesite()
    end
	if (para.logging && para.shifted) 		% log the shift if shifted (trustsite > 0)
		results.shift{para.loop} = para.shift;
	end
	if para.logging
		results.Vmat_svLog{para.loop} = cellfun(@(x) x(1,1), results.Vmat_sv(2:end));		% log the highest SV of Vmat
		results.EvaluesLog{para.loop} = results.Evalues;
	end
	if para.useFloShift3
		para.FloShift3minMaxSV = min(cellfun(@(x) x(1,1), results.Vmat_sv(2:end)));
	end
    fprintf('\nE = %.10g\t', results.E);							% print last energy Eigenvalue calculated
    results.Eerror(loop)=std(results.Evalues)/abs(mean(results.Evalues));
    fprintf('E_error =  %.16g\t',results.Eerror(loop));
    %Calculate the relative change of Vmat von Neumann entropy as one criteria for convergence check
    % TODO: Need to exclude Vmat_vNE == 0 as this gives either inf or NaN
    vNEdiff=(results.Vmat_vNE(2:end)-results.lastVmat_vNE(2:end))./results.Vmat_vNE(2:end);
    vNEdiff(results.Vmat_vNE(2:end)==0) = 0;        % if Vmat_vNE == 0 set diff to 0 to fix inf & NaN
    vNEdiff=abs(vNEdiff);
    fprintf('\nvNEdiff = \n'); disp(mat2str(vNEdiff));
    results.vNEdiff=vNEdiff;
    fprintf('para.shift = \n'); disp(mat2str(para.shift));
    results.lastVmat_vNE=results.Vmat_vNE;
    para=trustsite(para,results);
    fprintf('precise_sites k <= %d\t',para.precisesite);
    fprintf('trustable_sites k <= %d\t',para.trustsite(loop));

    if (para.dimlock==0 && para.trustsite(end)>3 && mod(loop,2)==0) || results.Eerror(end)<para.precision %&& para.trustsite(end)/para.precisesite>0.6 && sqrt(var(para.trustsite(end-2:end)))/mean(para.trustsite(end-2:end))<0.1 %The trust site has not been improved during the last 2 sweeps
         %%Expand or Truncate D and d_opt
         para.adjust=1;
%		 if results.Eerror(end)<para.precision						% PERHAPS GOOD STATEMENT!! TEST THIS! ADD 1 sweep after that before exit!
%			para.trustsite(end) = para.L;								% this ensures last adjustment of dimensions before finishing minimizeE.m
%		 end
         [op,para,results,mps,Vmat]=adjustdopt(op,para,results,mps,Vmat);
         [mps,Vmat,para, results] = rightnormA(mps,Vmat,para,results);				% changed to also update results
         %para.trustsite(loop)=0;
         fprintf('d_opt = ');
         disp(mat2str(para.d_opt));
         fprintf('para.D = ');
         disp(mat2str(para.D));
		 dispif('para.dk = ', para.useexpand);
         dispif(mat2str(para.dk),para.useexpand);
     else
        para.adjust=0;
        [mps,Vmat] = rightnormA(mps,Vmat,para,results);			%	sweeps r->l to right-normalize A matrices.
    end
    [op] = initstorage(mps, Vmat, op,para);							%	recreate operators for next optimization sweep

    fprintf('d_opt_change=%g\t',para.d_opt_change);
    fprintf('D_change=%g\t',para.D_change);
    if abs(para.d_opt_change)<para.minDimChange && abs(para.D_change)<para.minDimChange
      para.dimlock=1;
    end
    results.maxVmatsv=max(cellfun(@(x) x(end), results.Vmat_sv(2:para.trustsite(end))));fprintf('maxVmatsv=%g\t',results.maxVmatsv);
    results.maxAmatsv=max(cellfun(@(x) x(end), results.Amat_sv(2:para.trustsite(end))));fprintf('maxAmatsv=%g\n',results.maxAmatsv);
    %if para.trustsite(loop)>=para.L-2 && max(vNEdiff)<1e-3 && results.maxVmatsv<para.svmintol && results.maxAmatsv<para.svmintol
    if (results.Eerror(end)<para.precision && para.dimlock == 1) || (para.loop > 10 && prod(results.Eerror(end-10:end)<para.precision)) %&& para.trustsite(end)>para.L-5
		break;
    end

    %if para.trustsite(end)>para.L-10 && results.maxVmatsv<para.svmaxtol && results.maxAmatsv<para.svmaxtol
%	    break;
%    end
    %Save only every 10 sweeps to save time. (compared to save every step)
    if mod(loop,10)==0
    save(para.filename,'para','Vmat','mps','results','op');
    end
    loop=loop+1;
end
end
