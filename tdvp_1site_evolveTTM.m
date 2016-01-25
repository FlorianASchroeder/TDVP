function [mps,Vmat,para] = tdvp_1site_evolveTTM(mps,Vmat,para,results,op,outFile)
%% function tdvp_1site_evolveTTM(mps,Vmat,para,results,op,outFile)
% evolves the last site of the chain with TTM method.
% needs full tmps as input, needs storMPS = 1; read in MPS history from file!
%
% para.timeslice;		% current slice
Cl = cell(para.timeslice,1);				% stores bond-projectors for each past-MPS

% only load max. 100 MB at a time, calculate by taking current mps
vars = whos('mps','Vmat');
varSizeMB = sum([vars.bytes])/1024^2;
nPerBlock = floor(100/varSizeMB);			% n timeslices per block
currentN = 1;								% idx of start of not yet copied
while currentN < para.timeslice
	BlockN = min(nPerBlock, para.timeslice-currentN+1);
	tmps = outFile.tmps(currentN+(0:BlockN-1),:);
	tVmat = outFile.tVmat(currentN+(0:BlockN-1),:);
	
	% contract each slice in outFile with mps and Vmat
	for i = 1:BlockN
		Cl{currentN+i-1} = bondProject(i);		% i as index of tmps, tVmat
	end
	currentN = currentN + BlockN;
end


% only return on-site-MPS
mps = mps{end};
Vmat = Vmat{end};	

function Cl = bondProject(t)
%% function Cl = bondProject(timeslice)
% Projects the bond D(L-1) of mps(t*dt-1) onto D(L-1) of mps(tnow)
%
Cl = [];
% tmps  = outFile.tmps(t,:);			% read from file
% tVmat = outFile.tVmat(t,:);
for kk = 1:para.L-1
	Cl = updateCleft(Cl,mps{kk},Vmat{kk},[],tmps{t,kk},tVmat{t,kk});
end
% now Cl is Cl_(l',l) for site L; l': mps(tnow), l: mps(t*dt-1)
end

end

