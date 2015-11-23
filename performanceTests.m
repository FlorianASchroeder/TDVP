function out = performanceTests(variant, varargin)
%% function performanceTests(variant, varargin)
%	runs certain benchmarks which are useful to assess performance for TDVP-VMPS
out = [];
switch lower(variant)
	case 'matfile'
		testMatFile();
	case 'ncon'
		testNCON();
	case 'cellfun'
		testCellfun();
	case 'matcontract1'
		testGPUMatContract1();
	case 'tenscontract1'
		testGPUTensContract1();
	case 'tenscontract2'
		out = testGPUTensContract2();
end

end

function testMatFile()

%% Test performance of matfile command
timings = zeros(1,3);
for n = 10.^(1:3)
	system(['del "test.mat"'] );
	out = matfile('test.mat','Writable',true);
	start = tic;
	out.mps(n,n) = {[]};
	timings(log10(n)) = toc(start);
end
%%
	L = 100;
	N = [5,5,20];
	T = 1000;
    filespec = 'matfile_test.mat';
	%%
    mat = rand( N, 'double' );
	mps = cell(T,L);
	for i = 1:T
		for j = 1:L
			mps{i,j} = mat;
		end
	end
    save( filespec, 'mps','N','-v7.3' )
	clear mat mps
    %%
    obj = matfile( filespec ,'Writable', true );
    %%
    tic, mfm = obj.mps; toc
    tic, h5m = h5read( filespec, '/mps' ); toc
    %%
    dfm  = mfm-mat;
    d5m  = h5m-mat;
    max(abs(dfm(:)))
    max(abs(d5m(:)))
    %% column wise load
    tic, mfm = obj.mps( :, 1 ); toc
    tic, h5m = h5read( filespec, '/mps', [1,1], [N,1] ); toc
    %%
    dfm  = mfm-mat( :, 1 );
    d5m  = h5m-mat( :, 1 );
    max(abs(dfm(:)))
    max(abs(d5m(:)))
    %% row wise load
    tic, mfm = obj.mps( 2, : ); toc
    tic, h5m = h5read( filespec, '/mps', [2,1], [1,N] ); toc
	max(abs(mfm-h5m))
    %%
    dfm  = mfm-mat( 1, : );
    d5m  = h5m-mat( 1, : );
    max(abs(dfm(:)))
    max(abs(d5m(:)))
	%% writing test

	tic, obj.mat( 1, : ) = mfm; toc
%     tic, h5write( filespec, '/mat', h5m, [2,1], [1,N] ); toc

end

function testNCON()

%%  Benchmark updateCright vs ncon
tNCON = []; tUCR = []; step = 10;
for j = 4
	D = j*step; dOBB = 50; dk = 1000;
	mps = randn(D,D,dOBB);
	Vmat = randn(dk,dOBB);
	H = randn(dk,dk);
	C = randn(D,D);
	rep = 100;

	t1 = tic;
	for i = 1:rep
		cNCON =		ncon({mps,		conj(mps), Vmat,  conj(Vmat), H,	 C},...
			 {[-2,6,4], [-1,5,2],  [3,4], [1,2],	  [1,3], [5,6]});
	end
	tNCON(j) = toc(t1)

	t1 = tic;
	for i = 1:rep
		cUCR = updateCright(C, mps, Vmat, H, mps, Vmat);
	end

	tUCR(j) = toc(t1)
end

%% Plot results:
plot((1:5) .*step,[tUCR;tNCON]);

end

function testCellfun()
%% find(cellfun('isempty',mcOp)) --> fast enough!
rep = 1000;
mcOp = cell(10,1);
mcOp{4} = eye(10);
tic
for i = 1:rep
	ind = zeros(1,length(mcOp));
	ind = find(cellfun('isempty',mcOp));
end
toc

tic
for i = 1:rep
	ind = [];
	for j = 1:length(mcOp)
		if isempty(mcOp{j})
			ind = [ind;j];
		end
	end
end
toc

end

function testGPUMatContract1()
%% Test CPU vs GPU performance (matrix)
Nruns = 10;
dims = 10:20:5000;
tCPU = zeros(length(dims),Nruns);
tGPU = zeros(length(dims),Nruns);
gd = gpuDevice;

for ii = 1:length(dims)
	A = randn(dims(ii));
	B = randn(dims(ii));
	C = gpuArray(A);
	D = gpuArray(B);

	for jj = 1:Nruns
		t = tic;
		AB = A*B;
		tCPU(ii,jj) = toc(t);

		t = tic;
		CD = C*D;
		wait(gd);
		tGPU(ii,jj) = toc(t);
	end

	fprintf('D = %d, Average CPU time: %g, STD: %g\n', dims(ii), mean(tCPU(ii,:)), std(tCPU(ii,:)));
	fprintf('D = %d, Average GPU time: %g, STD: %g\n', dims(ii), mean(tGPU(ii,:)), std(tGPU(ii,:)));

end
end

function testGPUTensContract1()
%% Test CPU vs GPU performance (tensor)
Nruns = 10;
M = 6000;
N = 100;
P = 20;
K = 100;
dims = 10:20:K;
tCPU = zeros(length(dims),Nruns,4);
tGPU = zeros(length(dims),Nruns,4);
gd = gpuDevice;


for ii = 1:length(dims)
	A = rand(M,dims(ii),P);		% M K P
	B = rand(dims(ii),N);		% K N
	C = gpuArray(A);
	D = gpuArray(B);
	CD2 = gpuArray.zeros(M,N,P);

	for jj = 1:Nruns
		t = tic;
		AB = permute(reshape(reshape(permute(A,[1,3,2]),[],dims(ii))*B, M,P,N),[1,3,2]);
		tCPU(ii,jj,1) = toc(t);

		t = tic;
		ABcell = arrayfun(@(kk) A(:,:,kk) * B, 1:P,'UniformOutput',false);
		tCPU(ii,jj,2) = toc(t);

		t = tic;
		for kk = 1:P
			AB(:,:,kk) = A(:,:,kk) * B;
		end
		tCPU(ii,jj,3) = toc(t);

		% CONTRACTTENSORS
		t = tic;
		AB = permute(contracttensors(A,3,2,B,2,1),[1,3,2]);
		tCPU(ii,jj,4) = toc(t);

		% PAGEFUN
		t = tic;
		CD = pagefun(@mtimes, C,D);
		wait(gd);
		tGPU(ii,jj,1) = toc(t);

		% PERMUTE
		t = tic;
		CD = permute(reshape(reshape(permute(C,[1,3,2]),[],dims(ii))*D, M,P,N),[1,3,2]);
		wait(gd);
		tGPU(ii,jj,2) = toc(t);

		% FOR-LOOP
		t = tic;
		for kk = 1:P
			CD2(:,:,kk) = C(:,:,kk) * D;
		end
		wait(gd);
		tGPU(ii,jj,3) = toc(t);
		%parfor was pretty bad!

		% ARRAYFUN
		t = tic;
		CDcell = arrayfun(@(kk) C(:,:,kk) * D, 1:size(C,3),'UniformOutput',false);

		wait(gd);
		tGPU(ii,jj,4) = toc(t);

		% CONTRACTTENSORS
		t = tic;
		CD = permute(contracttensors(C,3,2,D,2,1),[1,3,2]);
		wait(gd);
		tGPU(ii,jj,5) = toc(t);
	end

	fprintf('D = %d, Average CPU time (permute): %g, STD: %g\n', dims(ii), mean(tCPU(ii,:,1)), std(tCPU(ii,:,1)));
	fprintf('D = %d, Average CPU time (arrayfu): %g, STD: %g\n', dims(ii), mean(tCPU(ii,:,2)), std(tCPU(ii,:,2)));
	fprintf('D = %d, Average CPU time (forloop): %g, STD: %g\n', dims(ii), mean(tCPU(ii,:,3)), std(tCPU(ii,:,3)));
	fprintf('D = %d, Average CPU time (contens): %g, STD: %g\n', dims(ii), mean(tCPU(ii,:,4)), std(tCPU(ii,:,4)));
	fprintf('D = %d, Average GPU time (pagefun): %g, STD: %g\n', dims(ii), mean(tGPU(ii,:,1)), std(tGPU(ii,:,1)));
	fprintf('D = %d, Average GPU time (permute): %g, STD: %g\n', dims(ii), mean(tGPU(ii,:,2)), std(tGPU(ii,:,2)));
	fprintf('D = %d, Average GPU time (forloop): %g, STD: %g\n', dims(ii), mean(tGPU(ii,:,3)), std(tGPU(ii,:,3)));
	fprintf('D = %d, Average GPU time (arrayfu): %g, STD: %g\n', dims(ii), mean(tGPU(ii,:,4)), std(tGPU(ii,:,4)));
	fprintf('D = %d, Average GPU time (contens): %g, STD: %g\n', dims(ii), mean(tGPU(ii,:,5)), std(tGPU(ii,:,5)));
	fprintf('\n');
end
	%%	Plot results
	f = figure(); f.Name = 'GPUbench'; ax=gca; hold all
	pl = {};
	pl{1} = plot(dims,mean(tCPU(:,:,1),2),'DisplayName','CPU, permute');
	pl{2} = plot(dims,mean(tCPU(:,:,2),2),'DisplayName','CPU, arrayfun');
	pl{3} = plot(dims,mean(tCPU(:,:,3),2),'DisplayName','CPU, for-loop');
	pl{4} = plot(dims,mean(tCPU(:,:,4),2),'DisplayName','CPU, con-tens');
	ax.ColorOrderIndex = 1;
	pl{5} = plot(dims,mean(tGPU(:,:,1),2),'-.','DisplayName','GPU, pagefun');
	pl{6} = plot(dims,mean(tGPU(:,:,2),2),'-.','DisplayName','GPU, permute');
	pl{7} = plot(dims,mean(tGPU(:,:,3),2),'-.','DisplayName','GPU, for-loop');
	pl{8} = plot(dims,mean(tGPU(:,:,4),2),'-.','DisplayName','GPU, arrayfun');
	pl{9} = plot(dims,mean(tGPU(:,:,5),2),'-.','DisplayName','GPU, con-tens');
	legend show
	xlabel('Array Dim');
	ylabel('Time in s');
	ax.Color = 'None';
	M*N*K*P

end

function out = testGPUTensContract2()
%% Test CPU vs GPU performance (tensor)
% test real StarmultA testcase
Nruns = 1000;
D = 20;		% max Dimension
dk = 4;		% local dimension
NC = 5;		% number of chains
d = [1, ones(1,NC)*D, dk];		% dimensions of A
dims = 1:D;
tCPU = nan(length(dims),Nruns);
tGPU = nan(length(dims),Nruns);
gd = gpuDevice;

for ii = 1:length(dims)
	d = [1, ones(1,NC)*dims(ii), dk];		% dimensions of A
	A = rand(d);
	op.h1jOBB = randn(dk);					% local single-site operator
	op.Hright = randn(dims(ii));			% right chain op
	op.Opright = randn(dims(ii));			% right chain op

% 	Agpu = gpuArray(A);

% 	CD2 = gpuArray.zeros(M,N,P);
	% CPU tensShape code
	jj = 0;
% 	while sum(tCPU(ii,:,1),'omitnan') < 20
	while (jj > 1 && (std(tCPU(ii,1:jj,1))/mean(tCPU(ii,1:jj,1)) > 0.1)) || jj < 5
		jj = jj + 1;
		t = tic;
		% 1. on-site H1
		w =	contracttensors(A, NC+2, NC+2, op.h1jOBB.', 2,1);

		for mc = 1:NC
			% Order for permute after contraction
% 			ord = [1:mc,NC+2,mc+1:NC+1];
			Atemp = tensShape(A, 'unfold', mc+1, d);		% chain index to front
			% 2. non-interacting Hlrstorage (Hright)
			OpTemp = op.Hright * Atemp;
			w = w + tensShape(OpTemp,'fold',mc+1, d);

			% 3. all interacting parts
			for mm = 1:2
				OpTemp = op.Opright * Atemp;					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
				OpTemp = tensShape(OpTemp,'fold',mc+1, d);
				w = w + contracttensors(OpTemp, NC+2, NC+2, op.h1jOBB.',2,1);
			end
		end
		tCPU(ii,jj,1) = toc(t);
		if jj == Nruns
			break;
		end
	end

	fprintf('D = %d, Average CPU time (tensShape): %g\n', dims(ii), mean(tCPU(ii,1:jj,1),'omitnan'));

	% GPU tensShape code (excluding transfer time)
	% copy all onto GPU
	Ag = gpuArray(A);
	op.h1jOBBg = gpuArray(op.h1jOBB);					% local single-site operator
	op.Hrightg = gpuArray(op.Hright);			% right chain op
	op.Oprightg = gpuArray(op.Opright);

	jj = 0;
% 	while sum(tCPU(ii,:,1),'omitnan') < 20
	while (jj > 1 && (std(tGPU(ii,1:jj,1))/mean(tGPU(ii,1:jj,1)) > 0.1)) || jj < 5
		jj = jj + 1;
		t = tic;
		% 1. on-site H1
		w =	contracttensors(Ag, NC+2, NC+2, op.h1jOBBg.', 2,1);

		for mc = 1:NC
			% Order for permute after contraction
% 			ord = [1:mc,NC+2,mc+1:NC+1];
			Atemp = tensShape(Ag, 'unfold', mc+1, d);		% chain index to front
			% 2. non-interacting Hlrstorage (Hright)
			OpTemp = op.Hrightg * Atemp;
			w = w + tensShape(OpTemp,'fold',mc+1, d);

			% 3. all interacting parts
			for mm = 1:2
				OpTemp = op.Oprightg * Atemp;					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
				OpTemp = tensShape(OpTemp,'fold',mc+1, d);
				w = w + contracttensors(OpTemp, NC+2, NC+2, op.h1jOBBg.',2,1);
			end
		end
		wait(gd);
		tGPU(ii,jj,1) = toc(t);
		if jj == Nruns
			break;
		end
	end

	fprintf('D = %d, Average GPU time (tShape-ex): %g\n', dims(ii), mean(tGPU(ii,1:jj,1),'omitnan'));

	% GPU tensShape code (including transfer time)
	jj = 0;
% 	while sum(tGPU(ii,:,2),'omitnan') < 20
	while (jj > 1 && (std(tGPU(ii,1:jj,2))/mean(tGPU(ii,1:jj,2)) > 0.1)) || jj < 5
		jj = jj + 1;
		t = tic;
		% copy all onto GPU

		Ag = gpuArray(A);
		op.h1jOBBg = gpuArray(op.h1jOBB);					% local single-site operator
		op.Hrightg = gpuArray(op.Hright);			% right chain op
		op.Oprightg = gpuArray(op.Opright);

		% 1. on-site H1
		w =	contracttensors(Ag, NC+2, NC+2, op.h1jOBBg.', 2,1);

		for mc = 1:NC
			% Order for permute after contraction
% 			ord = [1:mc,NC+2,mc+1:NC+1];
			Atemp = tensShape(Ag, 'unfold', mc+1, d);		% chain index to front
			% 2. non-interacting Hlrstorage (Hright)
			OpTemp = op.Hrightg * Atemp;
			w = w + tensShape(OpTemp,'fold',mc+1, d);

			% 3. all interacting parts
			for mm = 1:2
				OpTemp = op.Oprightg * Atemp;					% (m,2,1) should be the operator of site 2 in the effective left basis for system site 1
				OpTemp = tensShape(OpTemp,'fold',mc+1, d);
				w = w + contracttensors(OpTemp, NC+2, NC+2, op.h1jOBBg.',2,1);
			end
		end
		wait(gd);
		tGPU(ii,jj,2) = toc(t);
		if jj == Nruns
			break;
		end
		reset(gd);
	end

	fprintf('D = %d, Average GPU time (tShape-in): %g\n', dims(ii), mean(tGPU(ii,1:jj,2),'omitnan'));
	fprintf('\n');
end

	%%	Plot results
	f = figure(); f.Name = 'GPUbench'; ax=gca; hold all
	pl = {};
	pl{1} = plot(dims,mean(tCPU(:,:,1),2,'omitnan'),'DisplayName','CPU, tensShape');
% 	pl{2} = plot(dims,mean(tCPU(:,:,2),2),'DisplayName','CPU, arrayfun');
% 	pl{3} = plot(dims,mean(tCPU(:,:,3),2),'DisplayName','CPU, for-loop');
% 	pl{4} = plot(dims,mean(tCPU(:,:,4),2),'DisplayName','CPU, con-tens');
% 	ax.ColorOrderIndex = 1;
	pl{2} = plot(dims,mean(tGPU(:,:,1),2,'omitnan'),'-.','DisplayName','GPU, tensShape');
	pl{3} = plot(dims,mean(tGPU(:,:,2),2,'omitnan'),'-.','DisplayName','GPU, tensShape+transfer');
% 	pl{7} = plot(dims,mean(tGPU(:,:,3),2),'-.','DisplayName','GPU, for-loop');
% 	pl{8} = plot(dims,mean(tGPU(:,:,4),2),'-.','DisplayName','GPU, arrayfun');
% 	pl{9} = plot(dims,mean(tGPU(:,:,5),2),'-.','DisplayName','GPU, con-tens');
% 	legend show
	xlabel('Array Dim');
	ylabel('Time in s');
	ax.Color = 'None';
	set(gca,'YScale','log');

	out.tCPU = tCPU;
	out.tGPU = tGPU;
end


%% How to replace tensor contraction with bsxfun
