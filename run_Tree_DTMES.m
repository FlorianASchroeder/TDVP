function run_Tree_DTMES(tmax,dt,initState,CTShift)
%% function run_Tree_DTMES()
%
%	quick function to run TreeMPS DTMES calculation on noggin

fileName = VMPS_TreeMPS(...
			'model','DTMESclust4',...
			'evolveSysTrotter',0,...
			'L',18,...
			'tmax',tmax,...
			'deltaT',dt,...
			'extractObsInterval',dt,...
			'InitialState',initState,...
			'CTShift',CTShift,...
			'useDkExpand',1,...
			'Observables','.n.nd.na.x.xd.xa.dm.heff.',...
			'maxOBBDim',40,...
			'maxBondDim',[50,20],...
			'saveInterval',20);

load(fileName);

% If resuming:
%names = ls(sprintf('20160312-0432-5*-DPMES5-7C-Tree-v73TCMde9-L18CT%g%s/results-Till1500Step0.1v73-OBBmax60-Dmax20-expvCustom700-1core.mat',CTShift,initState));
%load(strtrim(names));

% Deserialize if needed
if ~isstruct(treeMPS)
%	hlp_deserialize_all();
	Vars = whos;				% deserialise if needed
	for ii = 1:size(Vars)
		if strcmp(Vars(ii).class,'uint8')
			eval(sprintf('%s = hlp_deserialize(%s);',Vars(ii).name,Vars(ii).name));
		end
	end
	resumed = 1;
	para.tdvp.resume = 1;		% switch on resume!
else
	resumed = 0;
end

if ~exist('tresults','var')
	tresults = [];
end

%% Copy to scratch for computation
if ~strcmp(computer,'PCWIN64')
	[~, name] = system('hostname');
	para.tdvp.hostname = strtrim(name);
	para.tdvp.scratchDir = '/scratch/fayns2/TDVPtemp/'; tempFold = fileparts(para.filename);
	currentDir = pwd;
	addpath(currentDir);
	if ~exist(para.tdvp.scratchDir,'dir')
		mkdir(para.tdvp.scratchDir);
	end
	mkdir([para.tdvp.scratchDir,tempFold]);
	copyfile(para.tdvp.filename,[para.tdvp.scratchDir,tempFold]);
	save(sprintf([para.tdvp.filename(1:end-4),'-incomplete-%s.mat'],para.tdvp.hostname),'para','results');
	cd(para.tdvp.scratchDir);
end

%% start time-evolution
para.tdvp.starttime = tic;
if para.useTreeMPS
	tdvp_1site_tree(treeMPS,para,results,tresults);
elseif para.useStarMPS
	tdvp_1site_star(mps,Vmat,para,results,op,tresults);
else
	tdvp_1site(mps,Vmat,para,results,op,tresults);
end

%% Copy Back to Home Folder
if ~strcmp(computer,'PCWIN64')
	copyfile(para.tdvp.filenameSmall,[currentDir,'/',para.tdvp.filenameSmall]);
	%delete([currentDir,'/',para.tdvp.filename(1:end-4),'-incomplete.mat']);
% 	sendmailCAM('fayns2@cam.ac.uk',...
%          'TDVP job completed',sprintf('The job \n %s\nHas successfully completed.',para.tdvp.filename));
	exit;
end

end
