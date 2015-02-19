function resume_TDVP(file)
% resumes a TDVP calculation
% Input: file = incomplete file in rscratch
maxNumCompThreads(1);						% safer in terms of results!
load(file);
[~, name] = system('hostname');
hostname = strtrim(name);
if ~strcmp(hostname, para.tdvp.hostname)
	error('You must continue job on the same PC!')
end
% remember current directory
currentDir = pwd;
addpath(currentDir);
% goto scratch and load file
cd(para.tdvp.scratchDir);
fprintf('Loading MPS...\n');
load(para.tdvp.filename);							% will contain the mps
fprintf('Loading parameters...\n');
load([para.tdvp.filename(1:end-4),'-small.mat']);	% all the rest
minSize = min(size(tmps,1),size(tVmat,1));	% perhaps -1 to be safer!
if size(tmps,1) ~= size(tVmat,1)
	% If there was a problem while saving state
	tmps    = tmps(1:minSize,:);
	tVmat   = tVmat(1:minSize,:);
end
if minSize >= tresults.lastIdx+1
	% If there is more MPS than calculated results
	calTimeFor = (tresults.lastIdx+1):minSize;
	fprintf('Calculating previously missing tresults...\n');
	tresults = calTimeObservables(tmps(calTimeFor,:),tVmat(calTimeFor,:),para,tresults);
end
fprintf('Saving MPS...\n');
save(para.tdvp.filename, 'tmps','tVmat','-v7.3');
fprintf('Saving parameters...\n');
save([para.tdvp.filename(1:end-4),'-small.mat'],'para','op','results','tresults','-v7.3');

mps  = tmps(end,:);
Vmat = tVmat(end,:);
clear('tmps','tVmat');
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