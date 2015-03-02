function resume_TDVP(file)
% resumes a TDVP calculation
% Input: file = incomplete file in rscratch
maxNumCompThreads(1);						% safer in terms of results!
load(file);
assert(isfield(para.tdvp,'filename'),'Need definition of para.tdvp.filename to save results!');

if ~strcmp(computer,'PCWIN64')
	[~, name] = system('hostname');
	hostname = strtrim(name);
	if ~strcmp(hostname, para.tdvp.hostname)
		error('You must continue job on the same PC!')
	end
	% remember current directory
	currentDir = pwd;
	addpath(currentDir);
	% goto scratch and load file
	if ~isfield(para.tdvp,'scratchDir')
		para.tdvp.scratchDir = '/scratch/fayns/TDVPtemp/';
	end
	cd(para.tdvp.scratchDir);
end

needSave = 0;
% Load files to calculate old tresults and initialize MPS and Vmat properly
fprintf('Loading MPS...\n');
outFile = matfile(para.tdvp.filename,'Writable',true);
if ~isprop(outFile,'tmps') || ~isprop(outFile,'tVmat')
	load(para.tdvp.filename);
	if ~(isfield(para.tdvp,'storeMPS') && para.tdvp.storeMPS == 0)
		error('There is no tMPS in the file');
	else
		clear('outFile');
	end
elseif size(fieldnames(outFile),1) > 4
	% more vars than tmps, tVmat, currentSize
	% -> save into new format
end
if exist('outFile','var')
	% Only pre-v41 files have not pre-allocated MPS-cells. Thus check sizes!
	[mtmps, ntmps] = size(outFile, 'tmps');
	[mtV,   ntV]   = size(outFile, 'tVmat');
	minSize = min(mtmps,mtV);						% perhaps -1 to be safer!
	if mtmps ~= mtV
		% If there was a problem while saving state
		outFile.tmps    = outFile.tmps(1:minSize,:);
		outFile.tVmat   = outFile.tVmat(1:minSize,:);
	end
end

fprintf('Loading parameters...\n');
try
	load([para.tdvp.filename(1:end-4),'-small.mat']);	% all the rest
catch err
	fprintf([getReport(err),'\n']);
	fprintf('\n');
end

if exist('outFile','var')
	% Pre-allocate MPS-cells for old files
	if mtmps ~= length(para.tdvp.t)
		fprintf('Allocating more space for MPS..');
		mtmps = length(para.tdvp.t);
		mtV = mtmps;
		outFile.tmps(mtmps,ntmps) = {[]};
		outFile.tVmat(mtV,ntV)  = {[]};
	end
	if ~isprop(outFile,'currentSize')
		outFile.currentSize = sum(~cellfun(@isempty,outFile.tmps(1:mtmps,1)));
	end
	minSize = outFile.currentSize;
	if minSize >= tresults.lastIdx+1
		% If there is more MPS than calculated results
		calTimeFor = (tresults.lastIdx+1):minSize;
		fprintf('Calculating previously missing tresults...\n');
		tmps = outFile.tmps(calTimeFor,:);				% only extract the needed parts of the MPS
		tVmat = outFile.tVmat(calTimeFor,:);
		tresults = calTimeObservables(tmps,tVmat,para,tresults);
		needSave = 1;
	end
	% take last slice to start calculations
	mps  = outFile.tmps(minSize,:);
	Vmat = outFile.tVmat(minSize,:);
	% clear('tmps','tVmat');
end

% fprintf('Saving MPS...\n');
% save(para.tdvp.filename, 'tmps','tVmat','-v7.3');
% 	expandBy = length(para.tdvp.t)-1 - size(results.tdvp.expError,1);
% 	results.tdvp.expError = [max(results.tdvp.expError,[],2); zeros(expandBy, 1)];	% take only the maximum error per sweep

if needSave
	fprintf('Saving parameters...\n');
	save([para.tdvp.filename(1:end-4),'-small.mat'],'para','op','results','tresults','-v7.3');
else
	fprintf('Resume TDVP without saving\n');
end

% resume the TDVP with the last matrices
para.tdvp.starttime = tic;
tdvp_1site(mps,Vmat,para,results,op);

% copy back to original dir
copyfile([para.tdvp.filename(1:end-4),'-small.mat'],[currentDir,'/',para.tdvp.filename(1:end-4),'-small.mat']);
%delete([currentDir,'/',para.tdvp.filename(1:end-4),'-incomplete.mat']);
%sendmailCAM('fayns2@cam.ac.uk',...
%	 'TDVP job completed',sprintf('The job \n %s\nHas successfully completed.',para.tdvp.filename));
return;

end