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

%% Test CPU vs GPU performance (tensor)
Nruns = 10;
dims = 10:20:2000;
tCPU = zeros(length(dims),Nruns);
tGPU = zeros(length(dims),Nruns);
gd = gpuDevice;

for ii = 1:length(dims)
	A = rand(20,dims(ii),50);
	B = rand(dims(ii),40);
	C = gpuArray(A);
	D = gpuArray(B);

	for jj = 1:Nruns
		t = tic;
		A = permute(A,[1,3,2]);
		A = reshape(A,[],dims(ii));
		AB = A*B;
		AB = reshape(AB,20,50,[]);
		AB = permute(AB,[1,3,2]);
		tCPU(ii,jj) = toc(t);

		t = tic;
		CD = pagefun(@mtimes, C,D);
		wait(gd);
		tGPU(ii,jj) = toc(t);
	end

	fprintf('D = %d, Average CPU time: %g, STD: %g\n', dims(ii), mean(tCPU(ii,:)), std(tCPU(ii,:)));
	fprintf('D = %d, Average GPU time: %g, STD: %g\n', dims(ii), mean(tGPU(ii,:)), std(tGPU(ii,:)));

end
	%%	Plot results
	f = figure(); f.Name = 'GPUbench'; ax=gca; hold all
	plot(dims,mean(tCPU,2));
	plot(dims,mean(tGPU,2));
	xlabel('Array Dim');
	ylabel('Time in s');
	ax.Color = 'None';