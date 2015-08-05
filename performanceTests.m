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

