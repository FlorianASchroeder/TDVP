%% Path where to save figures
path = pwd;
saveto = '..\Presentations\20140206 - Summary\';
%saveto = '..\..\Presentations\20140520 - MLSBM\20140525-CoupFunc\';
wantSave = 0;
%% Plot Vmat contributions for different Sites i (normalized)
i=9;
plot(real(Vmat{1,i}(:,:)))
title(['k = ',num2str(i),', max SV = ',num2str(results.Vmat_sv{1,i}(1,1))])
ylabel('Contribution to OBB')
xlabel('$d_k$')
if wantSave
    export_fig(sprintf('%sVmatNormalized%s-%u',saveto,para.filename(1:13),i),'-transparent','-png','-eps')
end

%% Plot Vmat contributions for different Sites i (SV-weighted)
i = 2;
plot(real(Vmat{1,i}*diag(results.Vmat_sv{1,i})))
title(['k = ',num2str(i),', max SV = ',num2str(results.Vmat_sv{1,i}(1,1))])
ylabel('Contribution to OBB')
xlabel('$d_k$')
%print(gcf, [saveto,'VmatScaled',num2str(i),'.eps'],'-deps')
if wantSave
    export_fig(sprintf('%sVmatScaled%s-%u',saveto,para.filename(1:13),i),'-transparent','-png','-eps')
end
%% Plot Sum over Vmat
i = 2;
a=sum(abs(real(Vmat{1,i}*diag(results.Vmat_sv{1,i}))),2);
plot(sum(abs(real(Vmat{1,i}*diag(results.Vmat_sv{1,i}))),2))
set(gca,'YScale','log')
title(['k = ',num2str(i),', max SV = ',num2str(results.Vmat_sv{1,i}(1,1))])
ylabel('Contribution to OBB')
xlabel('$d_k$')
%print(gcf, [saveto,'VmatScaled',num2str(i),'.eps'],'-deps')
if wantSave
    export_fig(sprintf('%sVmatScaled%s-%u',saveto,para.filename(1:13),i),'-transparent','-png','-eps')
end
%% Plot Sum over Vmat in 3D
plotMat = [];
rotate3d off
zeroVals = -50;             % value of zeros for padding and -inf replacement
for i = 2:length(Vmat)
    if size(Vmat{1,i},2) < size(results.Vmat_sv{1,i},1)
        a = log10(sum(abs(real(Vmat{1,i}*diag(results.Vmat_sv{1,i}(1:size(Vmat{1,i},2),:)))),2));
    else
        a = log10(sum(abs(real(Vmat{1,i}(:,1:size(results.Vmat_sv{1,i},1))*diag(results.Vmat_sv{1,i}))),2));
    end
    a(a==-inf)=zeroVals;
    %a = log10(sum(abs(real(Vmat{1,i}*diag(results.Vmat_sv{1,i}))),2));
    dim = max(length(a),size(plotMat,1));
    if length(a)< dim
        a = padarray(a,dim-length(a),zeroVals,'pre');
    elseif size(plotMat,1) < dim
        plotMat = padarray(plotMat,dim-size(plotMat,1),zeroVals,'pre');
    end
    plotMat = [plotMat,a];
end
surf(plotMat)
title(['k = ',num2str(i),', max SV = ',num2str(results.Vmat_sv{1,i}(1,1))])
ylabel('$d_k$')
xlabel('Site $k$')
set(gca,'View',[9.5 40]);
formatPlot(1)
rotate3d on

%% Plot Results
f1=figure(1);
for i = 2:1:size(results.D,2)
    if isempty(results.D{i})
        results.D{i} = results.D{i-1};
    end
    if isempty(results.d_opt{i})
        results.d_opt{i} = results.d_opt{i-1};
    end
    if isempty(results.dk{i})
        results.dk{i} = results.dk{i-1};
    end
end
for i = 2:1:size(results.shift,2)
    if isempty(results.shift{i})
        results.shift{i} = results.shift{i-1};
    end
end
subplot(2,2,1);
    plot(real(results.nx));
    title('$$<n_x(k)>$$');
subplot(2,2,2);
    plot(para.trustsite);
    title('Trustsite')
% 3D-version elucidating change:
if para.useshift
    subplot(2,2,3);
        surf(cell2mat(results.shift'))
        set(gca,'View',[-25 10]);
        shading interp
        title('Bosonic shift');
end
subplot(2,2,4);
    surf(cell2mat(results.d_opt'));
    shading interp
    set(gca,'View',[0 90]);
    title('OBB dim')
% 2D-version showing final
% subplot(2,2,3);
%     plot(results.shift{end});
%     title('Bosonic shift');
% subplot(2,2,4);
%     plot(results.d_opt{end});
%     title('OBB dim')
%     text(-80,-30,sprintf(para.filename(1:38)))
if wantSave
    export_fig(sprintf('%sResultsSummary%s',saveto,para.filename(1:13)),'-transparent','-png','-painters')
end

%% Plot d_opt, D adjustments
% fill cell arrays
for i = 2:1:size(results.D,2)
    if isempty(results.D{i})
        results.D{i} = results.D{i-1};
    end
    if isempty(results.d_opt{i})
        results.d_opt{i} = results.d_opt{i-1};
    end
    if isempty(results.dk{i})
        results.dk{i} = results.dk{i-1};
    end
end
f2 = figure(2);
subplot(1,3,1);
surf(cell2mat(results.D'))
set(gca,'View',[0 90]);
shading interp
title('Change in bond dimension D')
subplot(1,3,2);
surf(cell2mat(results.d_opt'))
set(gca,'View',[0 90]);
shading interp
title('Change in $d_{opt}$')
subplot(1,3,3);
surf(cell2mat(results.dk'))
set(gca,'View',[0 90]);
shading interp
title('Change in $d_{k}$')
rotate3d on
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure, to make caption readable
%export_fig(['png/',para.folder,'-D-dopt.png'],'-transparent',f2)
if wantSave
    export_fig(sprintf('%sChangeOfdoptDdk%s',saveto,para.filename(1:13)),'-transparent','-png')
end

%% Plot shift adjustment history of results.shift
saveto = '..\Presentations\20140206 - Summary\';
% fill cell array
for i = 2:1:size(results.shift,2)
    if isempty(results.shift{i})
        results.shift{i} = results.shift{i-1};
    end
end
f3 = figure(3);
%subplot(1,2,1);
%surf(cell2mat(results.D'))
%shading interp
%title('Change in bond dimension D')
%subplot(1,2,2);
surf(cell2mat(results.shift'))
set(gca,'View',[-25 10]);
shading interp
title('Change of shift $\delta_k$')
xlabel('site k')
ylabel('sweep')
zlabel('$\delta_k$')
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure, to make caption readable
if wantSave
    export_fig(sprintf('%sChangeOfShift%s',saveto,para.filename(1:13)),'-transparent','-png')
end

%% Try to analyse Differential Shift
diffShift = cell(1);
diffShift{1} = zeros(1,para.L);
diffShift{2} = zeros(1,para.L);
for i = 2:1:size(results.shift,2)
    if isempty(results.shift{i})
        results.shift{i} = results.shift{i-1};
    end
    if i~=2
        diffShift{i} = results.shift{i}-results.shift{i-1};
    end
end
f3 = figure(3);
%subplot(1,2,1);
%surf(cell2mat(results.D'))
%shading interp
%title('Change in bond dimension D')
%subplot(1,2,2);
surf(abs(cell2mat(diffShift')))
%set(gca,'View',[-25 10]);
set(gca,'View',[0 80]);
shading interp
title('Differential shift $\Delta (\delta_k)$')
xlabel('site k')
ylabel('sweep')
zlabel('$\delta_k$')
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure, to make caption readable
if wantSave
    export_fig(sprintf('%sDifferentialShift%s',saveto,para.filename(1:13)),'-transparent','-png')
end
%% plot results.Vmat_svLog
f4 = figure(4);
surf(cell2mat(results.Vmat_svLog'))
shading interp
title('Maximum SV of V')
xlabel('site k')
ylabel('sweep')
zlabel('SV')
set(gca,'View',[0 90]); %top
%set(gca,'View',[90 0]); %side
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure, to make caption readable
rotate3d on
if wantSave
    export_fig(sprintf('%sChangeOfVmatSV%s',saveto,para.filename(1:13)),'-transparent','-png')
end

%% plotting first rows of uneven cell arrays
cellfun(@(x) x(1,1), results.Vmat_sv(2:end))
%% try to calculate Shift analytically:
t = para.t; e = para.epsilon;
A = gallery('tridiag',t(2:end),e(1:end),t(2:end));    %creates tridiag for system A.(shiftVec) = (-t1*sigmaZ, 0,0,0,0,...)
B = zeros(para.L-1,1);
B(1) = -t(1)*1;%(results.spin.sz);  % -t1*sigmaZ
shift = A\B.*sqrt(2);
%% Plot my Shift versus VMPS Shift
colors={[0 0 1];[0 1 0];[1 0 0];[1 1 0];[1 0 1];[0 1 1];[0 0 0]};
hold all
pl(1) = plot([0;shift]);
set(pl(1), 'Marker', 'none','Color',colors{1}); % blue my
pl(2) = plot(para.shift);
set(pl(2), 'Marker', 'none','Color',colors{2}); % green VMPS
if wantSave
    export_fig(sprintf('%sAnalyticShift%s',saveto,para.filename(1:13)),'-transparent','-png','-eps')
end
%% Plot < n > of chain
figure(1)
pl(1) = plot(real(results.nx));
%set(gca,'YScale','log');
xlabel('Site k')
ylabel('$<n_{k,VMPS}>$')
formatPlot(1)
if wantSave
    export_fig(sprintf('%s%s-Occupation',saveto,para.filename(1:13)),'-transparent','-png','-pdf','-painters')
end
%% Plot Chain epsilon-t
figure(2)
title('Chain hopping and site energies');
subplot(1,2,1)
pl(1) = plot(para.t);
set(gca,'YScale','log');
xlabel('Site k')
ylabel('$t_k$')
subplot(1,2,2)
pl(2) = plot(para.epsilon);
set(gca,'YScale','log');
xlabel('Site k')
ylabel('$\epsilon_k$')

formatPlot(2)

%% Plot relative deviation from calculated < n >
%  Much more important as here the Wavefunction corrects also for different shift.
nx = [0 (shift.*shift./2)'];
figure(1)
plot((nx-real(results.nx))./nx,'LineStyle','none','Marker','*');
%set(gca,'YScale','log');
%title('Relative Deviation of $<n_k>$');
xlabel('Site k');
ylabel('$\frac{<n_k>-<n_{k,VMPS}>}{<n_k>}$')
formatPlot(1)
if wantSave
    export_fig(sprintf('%s%s-RelativeDeviationN',saveto,para.filename(1:13)),'-transparent','-png','-pdf','-painters')
end
%% Plot relative deviation of shift
relShift = ((shift-para.shift(2:end)')./shift);
figure(1)
pl(1) = plot([0; relShift],'LineStyle','none','Marker','*');
set(gca,'YScale','log');
%title('Relative deviation of shift from calculation')
xlabel('Site k')
ylabel('$\frac{\delta_k-\delta_{k,VMPS}}{\delta_k}$')
formatPlot(1)
if wantSave
    export_fig(sprintf('%sRelativeDeviationShift%s',saveto,para.filename(1:13)),'-transparent','-png','-pdf','-painters')
end
%% Plot ShiftUp vs ShiftDown
figure(5);
wantSave=0;
hold on
plot(results.bosonshift.x,'b')
plot(results.bosonshiftPerSpin.xUp,'r')
plot(results.bosonshiftPerSpin.xDown,'g')
title('Bosonic Shift vs Up or Down')
xlabel('site k')
ylabel('$\delta_k$')
legend('$x$','$x_\uparrow$','$x_\downarrow$');
ylim=get(gca,'YLim');
text(3,ylim(2)*2/3,...
    ['$x-(x_{\uparrow} +x_\downarrow) = $', sprintf('%.4g; \n', norm(results.bosonshift.x-results.bosonshiftPerSpin.xUp-results.bosonshiftPerSpin.xDown)),...
    '$|x|-|x_{\uparrow}|-|x_{\downarrow}| = $',sprintf('%.4g; ', norm(abs(results.bosonshift.x)-abs(results.bosonshiftPerSpin.xUp)-abs(results.bosonshiftPerSpin.xDown)))]);
formatPlot(5)
if wantSave
    export_fig(sprintf('%s%s-SBM-subOhmic-a0036622-shiftUpvsDown',saveto,para.filename(1:13)),'-transparent','-png','-painters')
end
%% Plot s-alpha relation
figure(5)
x = 0:0.01:1;
y = 0.005;
z = y.^(x./(1-x));
plot(x,z);

PlotData(:,PlotData(:,7)==0.001)

%% For PPC MLSBM:

%% Plot Energy convergence
figure(1);
plot(cell2mat(results.EvaluesLog)-min(cell2mat(results.EvaluesLog)));
disp(sprintf('%.15e',results.E))
set(gca,'YScale','log');
try
title(sprintf('$E_0 = %.10g, \\Lambda =  %.2g, z =  %.2g$',results.E, para.Lambda, para.z));catch end
xlabel('Site$\cdot$Loop');
ylabel('$E-E_0$');
formatPlot(1)
yLim = get(gca,'YLim');
for i = 1:para.loop
    line([para.L*i para.L*i],yLim,'LineWidth',1,'Color','black');
end
if wantSave
    export_fig(sprintf('%s%s-MLSBM-Econvergence-Lambda%.2gz%.2gp16',saveto,para.filename(1:13),para.Lambda,para.z),'-transparent','-png','-painters')
end

%% Get System-Bath coupling
figure(2);
plot(real(diag(op.h2term{1,1,1})./para.t(1)));
xlabel('Site k');
ylabel('$\hat\eta$');
formatPlot(2);

%% Plot System Wavefunction
figure(3);
plot(diag(calReducedDensity(mps,Vmat,para,1)))
if para.MLSB_staticDisorder
    hold all
    plot(para.MLSB_disDiag/500);        % /500 scales to see both
end
%%
%% Plot Flowdiagram
loop = para.loop;
plot(results.flowdiag{1,loop});

%% For TDVP analysis:

%% TDVP SBM: Plot evolution of the spin
figure(2);clf;
hold all
sphereon = true;
if sphereon
    sphere
    daspect([1 1 1])
    alpha(0.2)
end
col = parula(size(tresults.spin.sx,1));
scatter3(tresults.spin.sx,tresults.spin.sy,tresults.spin.sz,20,col,'filled');
% plot3(tresults.spin.sx,tresults.spin.sy,tresults.spin.sz);
set(gca,'xlim',[-1,1]);
set(gca,'ylim',[-1,1]);
set(gca,'zlim',[-1,1]);
rotate3d on

%% TDVP SBM: Plot Visibility / Coherence
figure(2);
% plot(para.tdvp.t, tresults.spin.visibility);
plot(para.tdvp.t(1:length(tresults.spin.sz)), tresults.spin.sz);
set(gca,'ylim',[-1,1]);
xlabel('t');
ylabel('$<s_z>$');

%% TDVP: Plot <n> environment
figure(3); clf;
n = size(tresults.nx,1);
surf(1:para.L,para.tdvp.t(1:n),real(tresults.nx))
xlabel('Site $k$');
ylabel('Time $t$');
zlabel('$<n_k>$');
% shading interp
rotate3d on
axis tight

%% TDVP: Plot temporal change in Vmat SV
figure(4); clf;
ax = axes('units','pixels');
pl = surf(cell2mat(results.tVmat_sv(1,:)));
set(gca,'View',[150 0]);
set(gca,'zscale','log');
set(gca,'zlimmode','manual');
set(gca,'zlim',[1e-15,1]);
xlabel('Site $k$');
ylabel('OBB Dimension');
% shading interp
sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',size(tmps,1),'Value',1,...
        'Position', [400 20 120 20],...
        'Callback', @(source,callbackdata) set(pl,'zdata',cell2mat(results.tVmat_sv(round(source.Value),:))));

%% TDVP: Plot Norm of Amat SV over time
figure(5);clf;
surf(1:(para.L-1),para.tdvp.t,cellfun(@norm, results.tdvp.Amat_sv))
xlabel('Site $k$');
ylabel('Time $t$');
axis tight
rotate3d on
%% TDVP: Plot Norm of Vmat SV over time
figure(5);clf;
surf(1:(para.L-1),para.tdvp.t,cellfun(@norm, results.tdvp.Vmat_sv(:,2:end)))
xlabel('Site $k$');
ylabel('Time $t$');
axis tight
rotate3d on

%% TDVP: Plot Sum over Vmat in 3D - analyse dk expand
figure(6);
rotate3d off
zeroVals = -50;             % value of zeros for padding and -inf replacement
plotMatTime = cell(size(tmps,1),1);
for j = 1:size(tmps,1)      % timeslices
    Vmatj = tVmat(j,:);
    Vmatj_sv = results.tdvp.Vmat_sv(j,:);
    plotMat = [];
for i = 2:length(Vmatj)
    if size(Vmatj{1,i},2) < size(Vmatj_sv{1,i},1)
        a = log10(sum(abs(real(Vmatj{1,i}*diag(Vmatj_sv{1,i}(1:size(Vmatj{1,i},2),:)))),2));
    else
        a = log10(sum(abs(real(Vmatj{1,i}(:,1:size(Vmatj_sv{1,i},1))*diag(Vmatj_sv{1,i}))),2));
    end
    a(a==-inf)=zeroVals;
    %a = log10(sum(abs(real(Vmat{1,i}*diag(results.Vmat_sv{1,i}))),2));
    dim = max(length(a),size(plotMat,1));
    if length(a)< dim
        a = padarray(a,dim-length(a),zeroVals,'pre');
    elseif size(plotMat,1) < dim
        plotMat = padarray(plotMat,dim-size(plotMat,1),zeroVals,'pre');
    end
    plotMat = [plotMat,a];
end
plotMatTime{j} = plotMat;
end
pl = surf(plotMatTime{1});
title('Temporal change of OBB')
ylabel('$d_k$')
xlabel('Site $k$')
set(gca,'View',[9.5 40]);
formatPlot(6)
axis tight
rotate3d on
sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',size(tmps,1),'Value',1,...
        'Position', [1100 20 120 20],...
        'Callback', @(source,callbackdata) set(pl,'zdata',plotMatTime{round(source.Value)}));

%% TDVP z-averaging in one file
% naming scheme to find files:
%   take series filename and replace z-value by *
folder = '20141025-1342-SpinBoson-alpha0.01delta0.1epsilon0dk20D5dopt5L49';
figure(7);clf;
% folder = '20141117-0531-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L84';
% folder = '20141117-0406-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L49';
filescheme = 'results-Till325Step4*-OBBExpand-noBondExpand*.mat';
files = dir(sprintf('%s/%s',folder,filescheme));
PlotData.spin.sz = [];
PlotData.z = [];
for k = 1:length(files)
    load([folder,'/',files(k).name]);
    PlotData.spin.sz(k,:) = tresults.spin.sz;
    PlotData.z(k) = para.z;
end
PlotData.t = para.tdvp.t;
plot(PlotData.t,[PlotData.spin.sz;mean(PlotData.spin.sz)]);
ylim([-1,1]);
legLabels = strsplit(sprintf('%.10g ',PlotData.z)); legLabels{end} = 'z-Ave';
legend(legLabels);
set(gca,'color','none');
xlabel('t');
ylabel('$<s_z>$');
%% Plot only z-averaged sz
plot(PlotData.t,mean(PlotData.spin.sz));
ylim([-1,1]);set(gca,'color','none');
xlabel('t');
ylabel('$<s_z>$');

%% TDVP (8) z-averaging create Orth 2010
% naming scheme to find files:
%   take series filename and replace z-value by *
figure(8);clf;
folder = {'20141025-1342-SpinBoson-alpha0.01delta0.1epsilon0dk20D5dopt5L49',...
          '20141114-2019-SpinBoson-alpha0.05delta0.1epsilon0dk20D5dopt5L49',...
          '20141117-0405-SpinBoson-alpha0.1delta0.1epsilon0dk20D5dopt5L49',...
          '20141117-0406-SpinBoson-alpha0.15delta0.1epsilon0dk20D5dopt5L49',...
          '20141117-0406-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L49'};
filescheme = 'results-Till325Step4*-OBBExpand-noBondExpand*.mat';

if ~exist('PlotData8','var')
    PlotData8.spin.meanSz = [];
    PlotData8.alpha = [];
end
if isempty(PlotData8.spin.meanSz)
    for k = 1:length(folder)
        PlotData8.spin.sz = [];
        files = dir(sprintf('%s/%s',folder{k},filescheme));
        for l = 1:length(files)
            load([folder{k},'/',files(l).name]);
            PlotData8.spin.sz(l,:) = tresults.spin.sz;
        end
        PlotData8.alpha(k) = para.alpha;
        PlotData8.spin.meanSz(k,:) = mean(PlotData8.spin.sz);
    end
end

PlotData8.t = para.tdvp.t;
plot(PlotData8.t,PlotData8.spin.meanSz);
ylim([-1,1]);
legLabels = strsplit(sprintf('%.10g ',PlotData8.alpha)); legLabels{end} = 'z-Ave';
legend(legLabels);
set(gca,'color','none');
xlabel('t');
ylabel('$<s_z>$');

%% TDVP (9) Orthogonal Polynomials L = 50: Orth 2010
figN = 10;
figure(figN); clf;
folder50 = {'20141114-1625-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50',...
            '20141114-1617-SpinBoson-OrthPol-alpha0.05delta0.1epsilon0dk20D5dopt5L50',...
            '20141117-0641-SpinBoson-OrthPol-alpha0.1delta0.1epsilon0dk20D5dopt5L50',...
            '20141117-0642-SpinBoson-OrthPol-alpha0.15delta0.1epsilon0dk20D5dopt5L50',...
            '20141116-0229-SpinBoson-OrthPol-alpha0.2delta0.1epsilon0dk20D5dopt5L50'};
if ~exist('PlotData9','var')
    PlotData9.sz = [];
    PlotData9.alpha = [];
end
if isempty(PlotData9.sz)
    for k = 1:length(folder50)
		% No OBB No Bond Expand:
%         filename = dir(sprintf('%s/results-Till325Step4-noOBBExpand-noBondExpand.mat',folder50{k}));
		% OBB and Bond Expand, max Bond 20;
		filename = dir(sprintf('%s/results-Till325Step4-OBBandBondExpand20*.mat',folder50{k}));
		try
	        load(sprintf('%s/%s',folder50{k},filename.name),'tresults','para');
		catch
			continue
		end
        PlotData9.sz(k,1:length(tresults.spin.sz)) = tresults.spin.sz;
        PlotData9.alpha(k) = para.alpha;
    end
    PlotData9.t = para.tdvp.t;
end
plot(PlotData9.t,PlotData9.sz);
ylim([-1,1]);
legLabels = strsplit(sprintf('%.10g ',PlotData9.alpha));
legend(legLabels(1:end-1));
set(gca,'color','none');
xlabel('t');
ylabel('$<s_z>$');
formatPlot(figN);

%% TDVP (10) Orthogonal Polynomials L = 200: Orth 2010
figN = 10;
figure(figN); clf;
folder50 = {'20141221-0148-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200',...
            '20141221-0148-SpinBoson-OrthPol-alpha0.05delta0.1epsilon0dk20D5dopt5L200',...
            '20141221-0151-SpinBoson-OrthPol-alpha0.1delta0.1epsilon0dk20D5dopt5L200',...
            '20141221-0151-SpinBoson-OrthPol-alpha0.15delta0.1epsilon0dk20D5dopt5L200',...
            '20141221-0153-SpinBoson-OrthPol-alpha0.2delta0.1epsilon0dk20D5dopt5L200'};
if ~exist('PlotData9','var')
    PlotData9.sz = [];
    PlotData9.alpha = [];
end
if isempty(PlotData9.sz)
    for k = 1:length(folder50)
		% OBB and Bond Expand, max Bond 20;
		filename = dir(sprintf('%s/results-Till325Step4-OBB*Expand20*.mat',folder50{k}))
		try
	        load(sprintf('%s/%s',folder50{k},filename.name),'tresults','para');
		catch
			continue
		end
        PlotData9.sz(k,1:length(tresults.spin.sz)) = tresults.spin.sz;
        PlotData9.alpha(k) = para.alpha;
    end
    PlotData9.t = para.tdvp.t;
end
plot(PlotData9.t,PlotData9.sz);
ylim([-1,1]);
legLabels = strsplit(sprintf('%.10g ',PlotData9.alpha));
legend(legLabels(1:end-1));
set(gca,'color','none');
xlabel('t');
ylabel('$<s_z>$');
formatPlot(figN);
%% TDVP (11) expvCustom Benchmarking
figure(9);clf;
hold all
scatter(sqrt(results.tdvp.expvTime(:,4)),results.tdvp.expvTime(:,1),'+');
scatter(sqrt(results.tdvp.expvTime(:,4)),results.tdvp.expvTime(:,2)+results.tdvp.expvTime(:,3),'*');
% scatter(sqrt(results.tdvp.expvTime(:,4)),results.tdvp.expvTime(:,2),'*');
% scatter(sqrt(results.tdvp.expvTime(:,4)),results.tdvp.expvTime(:,3));
set(gca,'yscale','log')
set(gca,'xscale','log')
legend('Custom Krylov e^{At}v','Expokit')
xlabel('Matrix dimension n')
ylabel('Time/s')
formatPlot(9)
