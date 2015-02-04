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
figure(1); hold on;
pl(1) = plot(real(results.nx));
%set(gca,'YScale','log');
xlabel('Site k')
ylabel('$<n_{k,VMPS}>$')
set(gca,'yscale','log')
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
%     line([para.L*i para.L*i],yLim,'LineWidth',1,'Color','black');
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
	set(get(gca,'children'),'linestyle',':')
end
col = parula(size(tresults.spin.sx,1));
scatter3(tresults.spin.sx,tresults.spin.sy,tresults.spin.sz,20,col,'filled');
plot3(tresults.spin.sx,tresults.spin.sy,tresults.spin.sz);
set(gca,'xlim',[-1,1]);
set(gca,'ylim',[-1,1]);
set(gca,'zlim',[-1,1]);
set(gca,'view',[-29,16]);
rotate3d on

%% TDVP SBM: Plot Visibility / Coherence
figure(2); hold all;
% plot(para.tdvp.t, tresults.spin.visibility);
plot(para.tdvp.t(1:length(tresults.spin.sz)), tresults.spin.sz);
set(gca,'ylim',[-1,1]);
xlabel('t');
ylabel('$<s_z>$');

%% TDVP: Plot <n> environment
figure(3); clf;
n = size(tresults.nx,1);
surf(1:para.L,para.tdvp.t(1:n),real(tresults.nx))
% surf(1:para.L,para.tdvp.t(1:n),abs(real(tresults.nx)-ones(size(tresults.nx,1),1)*real(tresults.nx(1,:))))
xlabel('Site $k$');
ylabel('Time $t$');
zlabel('$<n_k>$');
set(gca,'zscale','log');
set(gca,'View',[0 42]);
% shading interp
rotate3d on
axis tight

%% TDVP: Draw <n> propagation
figure(3); clf;
n = size(tresults.nx,1); nCont = 30;
contourf(1:para.L,para.tdvp.t(1:n),real(tresults.nx),nCont)

coneParam = zeros(2,n);
hold on
for m=1:8
	[coneParam(1,m),coneParam(2,m)] = ginput(1);
	plot(coneParam(1,m),coneParam(2,m),'ko','markerSize',10,'lineWidth',1,'markerFaceColor','w');
end
%%
pf = polyfit(coneParam(1,:),coneParam(2,:),3);
coneFunct = polyval(pf,para.tdvp.t(1:n));
plot(coneFunct,para.tdvp.t(1:n),'w','linewidth',4)
plot(coneFunct,para.tdvp.t(1:n),'k','linewidth',2)
hold off
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
folder = '20150131-1447-SpinBoson-LogZ-v37TCM67-alpha0.2delta0.1epsilon0dk20D5dopt5L30-artificial';
figure(7);clf;
% folder = '20141117-0531-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L84';
% folder = '20141117-0406-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L49';
filescheme = 'results-Till325Step4*-OBBExpand-noBondExpand*-small.mat';
files = ls(sprintf('%s/%s',folder,filescheme));
PlotData.spin.sz = [];PlotData.spin.sx = []; PlotData.spin.sy = []; PlotData.nx = [];
PlotData.z = [];
for k = 1:size(files,1)
    load([folder,'/',files(k,:)]);
	PlotData.spin.sx(k,:) = tresults.spin.sx;
	PlotData.spin.sy(k,:) = tresults.spin.sy;
    PlotData.spin.sz(k,:) = tresults.spin.sz;
	PlotData.nx(:,:,k)    = tresults.nx;
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
para.tdvp.filename = [para.tdvp.filename(1:end-5),'Avg.mat'];
tresults.spin.sx = mean(PlotData.spin.sx,1);
tresults.spin.sy = mean(PlotData.spin.sy,1);
tresults.spin.sz = mean(PlotData.spin.sz,1);
tresults.nx = mean(PlotData.nx,3);
if ~exist(para.tdvp.filename,'file')
	save(para.tdvp.filename, 'para','results','tresults');
end

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

%% TDVP (12) prepare artificial GroundState
% clear Workspace, load old result.mat and execute. Saves into new folder
% for fresh start! Initializes Spin in +Sz
para.folder = [para.folder,'-artificial'];
para.filename = sprintf('%s/results.mat',para.folder);
mps{1} = reshape([1,zeros(1,numel(mps{1})-1)],[1,para.D(1),para.d_opt(1)]);
Vmat{1} = eye(para.dk(1));
for j=2:para.L-1
	mps{j} = reshape([1, zeros(1,numel(mps{j})-1)],para.D(j-1),para.D(j),para.d_opt(j));
	Vmat{j} = [zeros(para.dk(j)-para.d_opt(j),para.d_opt(j));...
		fliplr(eye(para.d_opt(j)))];
end
mps{para.L} = reshape([1, zeros(1,numel(mps{para.L})-1)],para.D(para.L-1),1,para.d_opt(para.L));
Vmat{para.L} = [zeros(para.dk(para.L)-para.d_opt(para.L),para.d_opt(para.L));...
		fliplr(eye(para.d_opt(para.L)))];
% Vmat{para.L}(para.dk(para.L):-1:para.dk(para.L)-para.d_opt(para.L)+1,1:para.d_opt(para.L)) = -eye(para.d_opt(para.L));

% results.nx         = getObservable({'occupation'},mps,Vmat,para);
% results.bosonshift = getObservable({'shift'},mps,Vmat,para);

if strcmp(para.model,'SpinBoson')
%     results.spin   = getObservable({'spin'},mps,Vmat,para);
end
% mkdir(para.folder);
% save(para.filename,'mps','Vmat','para','results','op');

%% TDVP SBM multi load files: OrthPol rev22 threading vs perfect L=200
res = {};
res{1,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150115-1539-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand15-expvCustom800-1pass-small.mat');
res{1,2} = 'v22 multi-core, 1pass, rescaling=0';
res{2,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150115-1539-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand15-expvCustom800-1pass-small.mat');
res{2,2} = 'v22 multi-core, restarted, rescaling=0';
res{3,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150115-1504-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand15-expvCustom800-1pass-small.mat');
res{3,2} = 'v22 single-core, 1pass, rescaling=1';
res{4,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150115-1504-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand15-expvCustom800-small.mat');
res{4,2} = 'v22 single-core, restarted, rescaling=1';
res{5,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150115-1510-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand15-expvCustom800-1pass-small.mat');
res{5,2} = 'v22 single-core, 1pass, rescaling=0';
res{6,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150115-1510-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand15-expvCustom800-small.mat');
res{6,2} = 'v22 single-core, restarted, rescaling=0';
res{7,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150115-1539-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v22-OBBExpand-BondExpand25-expvCustom800-small.mat');
res{7,2} = 'v22 multi-core, restarted, rescaling=0, Bond25';
res{8,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpandBondExpand-small.mat');
res{8,2} = 'perfect';
res{9,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpand-BondExpand15-expvCustom800-small.mat');
res{9,2} = 'v22,expvCustom800,Bond15';
res{10,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpand-BondExpand15-small.mat');
res{10,2}= 'v22,noExvpVCustom,Bond15';
res{11,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1237-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{11,2}= 'v35, rescaling = 1, Lappy';
res{12,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1326-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{12,2}= 'v35, rescaling = 0, Lappy';
res{13,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1431-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{13,2}= 'v35, rescaling = 1, TCM';
res{14,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1432-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{14,2}= 'v35, rescaling = 0, TCM';
res{15,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1509-SpinBoson-OrthPol-exp.v20-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{15,2}= 'v20, rescaling = 1, TCM';
res{16,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1510-SpinBoson-OrthPol-exp.v20-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{16,2}= 'v20, rescaling = 0, TCM';


% cell2mat(cellfun(@(x) [x.results.time,x.results.tdvp.time], res(:,1), 'UniformOutput', false))

%% TDVP SBM multi load files: OrthPol L=200 rev20 vs rev35 GS
res{1,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpandBondExpand-small.mat');
res{1,2} = 'perfect v19';
res{2,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1509-SpinBoson-OrthPol-exp.v20-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{2,2}= 'v20, rescaling = 1, TCM';
res{3,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1510-SpinBoson-OrthPol-exp.v20-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{3,2}= 'v20, rescaling = 0, TCM';
res{4,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1602-SpinBoson-OrthPol-v26TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc1\results.mat');
res{4,2}= 'v26, rescaling = 1, TCM';
res{5,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1602-SpinBoson-OrthPol-v26TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc0\results.mat');
res{5,2}= 'v26, rescaling = 0, TCM';
res{6,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1554-SpinBoson-OrthPol-v27TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc1\results.mat');
res{6,2}= 'v27, rescaling = 1, TCM';
res{7,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1555-SpinBoson-OrthPol-v27TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc0\results.mat');
res{7,2}= 'v27, rescaling = 0, TCM';
res{8,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1607-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc1\results.mat');
res{8,2}= 'v28, rescaling = 1, TCM';
res{9,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1609-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc0\results.mat');
res{9,2}= 'v28, rescaling = 0, TCM';
res{10,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1614-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc1\results.mat');
res{10,2}= 'v30, rescaling = 1, TCM';
res{11,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1612-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc0\results.mat');
res{11,2}= 'v30, rescaling = 0, TCM';
res{12,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1237-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{12,2}= 'v35, rescaling = 1, Lappy';
res{13,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1326-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{13,2}= 'v35, rescaling = 0, Lappy';
res{14,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1431-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{14,2}= 'v35, rescaling = 1, TCM';
res{15,1}= load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1432-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results.mat');
res{15,2}= 'v35, rescaling = 0, TCM';

%% TDVP SBM multi load files: OrthPol compare L = 50 rev20, rev22 vs
res = {};
res{1,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20141114-1625-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4-noOBBExpand-noBondExpand.mat');
res{1,2} = '14/11/2014, noExpands, Exp.v19?';
res{2,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150120-1952-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results.mat');
res{2,2} = 'GS Exp.v25, rescaling = 1';
res{3,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP2\20150120-2047-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results.mat');
res{3,2} = 'GS Exp.v23, rescaling = 1';
res{4,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP2\20150120-2034-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results.mat');
res{4,2} = 'GS Exp.v22, rescaling = 1';
res{5,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP2\20150120-2059-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results.mat');
res{5,2} = 'GS Exp.v21, rescaling = 1';
res{6,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP2\20150120-2118-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results.mat');
res{6,2} = 'GS Exp.v20, rescaling = 1';
res{7,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150120-2255-SpinBoson-OrthPol-rev27-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{7,2} = 'GS to22.v27, rescaling = 1';
res{8,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150120-2310-SpinBoson-OrthPol-rev28-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{8,2} = 'GS to22.v28, rescaling = 1';
res{9,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150120-2321-SpinBoson-OrthPol-rev29-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{9,2} = 'GS to22.v29, rescaling = 1';
res{10,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150120-2346-SpinBoson-OrthPol-rev29-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{10,2} = 'GS to22.v29, rescaling = 1, higher precision';
res{11,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150121-0138-SpinBoson-OrthPol-rev29-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{11,2} = 'GS to22.v29, rescaling = 0';
res{12,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150121-0156-SpinBoson-OrthPol-rev28-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{12,2} = 'GS to22.v28, rescaling = 0';
res{13,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150121-1602-SpinBoson-OrthPol-rev30-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{13,2} = 'GS to22.v30, rescaling = 1';
res{14,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150121-1634-SpinBoson-OrthPol-rev30-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{14,2} = 'GS to22.v30, rescaling = 0';
res{15,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150121-1817-SpinBoson-OrthPol-rev31-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{15,2} = 'GS to22.v31, rescaling = 1';
res{16,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150121-1821-SpinBoson-OrthPol-rev31-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{16,2} = 'GS to22.v31, rescaling = 0';
res{17,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150121-2152-SpinBoson-OrthPol-rev32-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{17,2} = 'GS to22.v32, rescaling = 1';
res{18,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150121-2203-SpinBoson-OrthPol-rev32-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{18,2} = 'GS to22.v32, rescaling = 1, no braces';
res{19,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150122-0006-SpinBoson-OrthPol-rev33-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{19,2} = 'GS to22.v33, rescaling = 1, HmultA sum';
res{20,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150122-0018-SpinBoson-OrthPol-rev34-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{20,2} = 'GS to22.v34, rescaling = 1';
res{21,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1444-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{21,2} = 'GS Exp.v35, rescaling = 1';
res{22,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1447-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{22,2} = 'GS Exp.v35, rescaling = 0';
res{23,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1436-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{23,2} = 'GS Exp.v35, rescaling = 1, TCM';
res{24,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1438-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{24,2} = 'GS Exp.v35, rescaling = 0, TCM';
res{25,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1507-SpinBoson-OrthPol-exp.v20-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{25,2} = 'GS Exp.v20, rescaling = 1, TCM';
res{26,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150123-1508-SpinBoson-OrthPol-exp.v20-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{26,2} = 'GS Exp.v20, rescaling = 0, TCM';
res{27,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1600-SpinBoson-OrthPol-v26TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc1/results.mat');
res{27,2} = 'GS Exp.v26, rescaling = 1, TCM';
res{28,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1600-SpinBoson-OrthPol-v26TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc0/results.mat');
res{28,2} = 'GS Exp.v26, rescaling = 0, TCM';
res{29,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1552-SpinBoson-OrthPol-v27TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc1/results.mat');
res{29,2} = 'GS Exp.v27, rescaling = 1, TCM';
res{30,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1553-SpinBoson-OrthPol-v27TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc0/results.mat');
res{30,2} = 'GS Exp.v27, rescaling = 0, TCM';
res{31,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1605-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc1/results.mat');
res{31,2} = 'GS Exp.v28, rescaling = 1, TCM';
res{32,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1608-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc0/results.mat');
res{32,2} = 'GS Exp.v28, rescaling = 0, TCM';
res{33,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1612-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc1/results.mat');
res{33,2} = 'GS Exp.v30, rescaling = 1, TCM';
res{34,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1611-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc0/results.mat');
res{34,2} = 'GS Exp.v30, rescaling = 0, TCM';
res{35,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150126-1511-SpinBoson-OrthPol-exp.v35-alpha0.01delta0.1epsilon0dk20D5dopt5L50/results.mat');
res{35,2} = 'GS Exp.v35, rescaling = 1, TCM,prec';

%% TDVP SBM multi load files: OrthPol L= 50 rev28, rev30 architecture sweep
res{1,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1605-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc1/results.mat');
res{1,2} = 'GS Exp.v28, rescaling = 1, Haswell';
res{2,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1608-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc0/results.mat');
res{2,2} = 'GS Exp.v28, rescaling = 0, Haswell';
res{3,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1612-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc1/results.mat');
res{3,2} = 'GS Exp.v30, rescaling = 1, Haswell';
res{4,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1611-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L50resc0/results.mat');
res{4,2} = 'GS Exp.v30, rescaling = 0, Haswell';
res{5,1} = load([ls('201501*OrthPol*v28TCM*Sandy*L50resc1'),'/results.mat']);
res{5,2} = 'GS Exp.v28, rescaling = 1, Sandy';
res{6,1} = load([ls('201501*OrthPol*v28TCM*Sandy*L50resc0'),'/results.mat']);
res{6,2} = 'GS Exp.v28, rescaling = 0, Sandy';
res{7,1} = load([ls('201501*OrthPol*v30TCM*Sandy*L50resc1'),'/results.mat']);
res{7,2} = 'GS Exp.v30, rescaling = 1, Sandy';
res{8,1} = load([ls('201501*OrthPol*v30TCM*Sandy*L50resc0'),'/results.mat']);
res{8,2} = 'GS Exp.v30, rescaling = 0, Sandy';
res{9,1} = load([ls('201501*OrthPol*v28TCM*Nehalem*L50resc1'),'/results.mat']);
res{9,2} = 'GS Exp.v28, rescaling = 1, Nehalem';
res{10,1} = load([ls('201501*OrthPol*v28TCM*Nehalem*L50resc0'),'/results.mat']);
res{10,2} = 'GS Exp.v28, rescaling = 0, Nehalem';
res{11,1} = load([ls('201501*OrthPol*v30TCM*Nehalem*L50resc1'),'/results.mat']);
res{11,2} = 'GS Exp.v30, rescaling = 1, Nehalem';
res{12,1} = load([ls('201501*OrthPol*v30TCM*Nehalem*L50resc0'),'/results.mat']);
res{12,2} = 'GS Exp.v30, rescaling = 0, Nehalem';
res{13,1} = load([ls('201501*OrthPol*v28TCM*Core2*L50resc1'),'/results.mat']);
res{13,2} = 'GS Exp.v28, rescaling = 1, Core2';
res{14,1} = load([ls('201501*OrthPol*v28TCM*Core2*L50resc0'),'/results.mat']);
res{14,2} = 'GS Exp.v28, rescaling = 0, Core2';
res{15,1} = load([ls('201501*OrthPol*v30TCM*Core2*L50resc1'),'/results.mat']);
res{15,2} = 'GS Exp.v30, rescaling = 1, Core2';
res{16,1} = load([ls('201501*OrthPol*v30TCM*Core2*L50resc0'),'/results.mat']);
res{16,2} = 'GS Exp.v30, rescaling = 0, Core2';
%% TDVP SBM multi load files: OrthPol L=200 rev28, rev30 architecture sweep
res{1,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1607-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc1/results.mat');
res{1,2} = 'GS Exp.v28, rescaling = 1, Haswell';
res{2,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1609-SpinBoson-OrthPol-v28TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc0/results.mat');
res{2,2} = 'GS Exp.v28, rescaling = 0, Haswell';
res{3,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1614-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc1/results.mat');
res{3,2} = 'GS Exp.v30, rescaling = 1, Haswell';
res{4,1} = load('E:\Documents\Uni\PhD\Theory\schroederflorian-vmps-tdvp\TDVP\20150124-1612-SpinBoson-OrthPol-v30TCM-alpha0.01delta0.1epsilon0dk20D5dopt5L200resc0/results.mat');
res{4,2} = 'GS Exp.v30, rescaling = 0, Haswell';
res{5,1} = load([ls('201501*OrthPol*v28TCM*Sandy*L200resc1'),'/results.mat']);
res{5,2} = 'GS Exp.v28, rescaling = 1, Sandy';
res{6,1} = load([ls('201501*OrthPol*v28TCM*Sandy*L200resc0'),'/results.mat']);
res{6,2} = 'GS Exp.v28, rescaling = 0, Sandy';
res{7,1} = load([ls('201501*OrthPol*v30TCM*Sandy*L200resc1'),'/results.mat']);
res{7,2} = 'GS Exp.v30, rescaling = 1, Sandy';
res{8,1} = load([ls('201501*OrthPol*v30TCM*Sandy*L200resc0'),'/results.mat']);
res{8,2} = 'GS Exp.v30, rescaling = 0, Sandy';
res{9,1} = load([ls('201501*OrthPol*v28TCM*Nehalem*L200resc1'),'/results.mat']);
res{9,2} = 'GS Exp.v28, rescaling = 1, Nehalem';
res{10,1} = load([ls('201501*OrthPol*v28TCM*Nehalem*L200resc0'),'/results.mat']);
res{10,2} = 'GS Exp.v28, rescaling = 0, Nehalem';
res{11,1} = load([ls('201501*OrthPol*v30TCM*Nehalem*L200resc1'),'/results.mat']);
res{11,2} = 'GS Exp.v30, rescaling = 1, Nehalem';
res{12,1} = load([ls('201501*OrthPol*v30TCM*Nehalem*L200resc0'),'/results.mat']);
res{12,2} = 'GS Exp.v30, rescaling = 0, Nehalem';
res{13,1} = load([ls('201501*OrthPol*v28TCM*Core2*L200resc1'),'/results.mat']);
res{13,2} = 'GS Exp.v28, rescaling = 1, Core2';
res{14,1} = load([ls('201501*OrthPol*v28TCM*Core2*L200resc0'),'/results.mat']);
res{14,2} = 'GS Exp.v28, rescaling = 0, Core2';
res{15,1} = load([ls('201501*OrthPol*v30TCM*Core2*L200resc1'),'/results.mat']);
res{15,2} = 'GS Exp.v30, rescaling = 1, Core2';
res{16,1} = load([ls('201501*OrthPol*v30TCM*Core2*L200resc0'),'/results.mat']);
res{16,2} = 'GS Exp.v30, rescaling = 0, Core2';

%% TDVP SBM multi load files: OrthPol/LogZ v37: proper L=50, L=200 Sweep!
% OrthPol-TDVP-Benchmark-Sz-L50-v37:				1,3:6
% OrthPol-TDVP-Benchmark-Sz-L200-v37:				1,7:9
% OrthPol-TDVP-Benchmark-Sz-L50-artificial-v37:		1,10:14
% OrthPol-TDVP-Benchmark-Sz-L200-artificial-v37:	1,15:18
% Legend font size: 16 or 18
res{1,1} = load('20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpandBondExpand-small.mat');
res{1,2} = 'v19, Perfect';
res{2,1} = load('20150126-1800-SpinBoson-LogZ-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L49/results-Till325Step4v37-OBBExpand-noBondExpand-1core-zAvg.mat');
res{2,2} = '\Lambda=2, z_{Avg}=0.2, OBB no Bond';
res{3,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4v37-noOBBExpand-noBondExpand-1core-small.mat');
res{3,2} = 'OrthPol L=50, no Expand';
res{4,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4v37-OBBExpand-noBondExpand-1core-small.mat');
res{4,2} = 'OrthPol L=50, OBB no Bond';
res{5,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4v37-noOBBExpand-BondExpand15-1core-small.mat');
res{5,2} = 'OrthPol L=50, Bond no OBB';
res{6,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom800-1core-small.mat');
res{6,2} = 'OrthPol L=50, All Expand';
res{7,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v37-noOBBExpand-noBondExpand-1core-small.mat');
res{7,2} = 'OrthPol L=200, no Expand';
res{8,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v37-OBBExpand-noBondExpand-1core-small.mat');
res{8,2} = 'OrthPol L=200, OBB no Bond';
res{9,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom800-1core-small.mat');
res{9,2} = 'OrthPol L=200, All Expand';
res{10,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-noOBBExpand-noBondExpand-1core-small.mat');
res{10,2} = 'OrthPol L=50, Art, no Expand';
res{11,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-OBBExpand-noBondExpand-1core-small.mat');
res{11,2} = 'OrthPol L=50, Art, OBB no Bond';
res{12,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-noOBBExpand-BondExpand15-1core-small.mat');
res{12,2} = 'OrthPol L=50, Art, Bond no OBB';
res{13,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom800-1core-small.mat');
res{13,2} = 'OrthPol L=50, Art, All, expvCustom800';
res{14,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom1-1core-small.mat');
res{14,2} = 'OrthPol L=50, Art, All, expvCustom1';
res{15,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial\results-Till325Step4v37-noOBBExpand-noBondExpand-1core-small.mat');
res{15,2} = 'OrthPol L=200,Art,no Expand';
res{16,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial\results-Till325Step4v37-OBBExpand-noBondExpand-1core-small.mat');
res{16,2} = 'OrthPol L=200,Art, OBB no Bond';
res{17,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial\results-Till325Step4v37-noOBBExpand-BondExpand15-1core-small.mat');
res{17,2} = 'OrthPol L=200,Art, Bond no OBB';
res{18,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom800-1core-small.mat');
res{18,2} = 'OrthPol L=200,Art, All';
res{19,1} = load('20150126-2247-SpinBoson-LogZ-v37TCM66-alpha0.2delta0.1epsilon0dk20D5dopt5L49/results-Till325Step4v37-OBBExpand-noBondExpand-1core-zAvg.mat');
res{19,2} = '\Lambda=2, z_{Avg}=0.2, OBB no Bond';
res{20,1} = load('20150126-2247-SpinBoson-LogZ-v37TCM66-alpha0.2delta0.1epsilon0dk20D5dopt5L49-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-1core-zAvg.mat');
res{20,2} = '\Lambda=2, z_{Avg}=0.2, OBB no Bond, Art';

%% TDVP SBM multi load files: Orth2010, OrthPol, artificial, L=50, L=200
res{ 1,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 1,2} = '\alpha = 0.01';
res{ 2,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.05delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 2,2} = '\alpha = 0.05';
res{ 3,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.1delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 3,2} = '\alpha = 0.1';
res{ 4,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.15delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 4,2} = '\alpha = 0.15';
res{ 5,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.2delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 5,2} = '\alpha = 0.2';
res{ 6,1} = load('20150130-1531-SpinBoson-OrthPol-v37TCM33-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 6,2} = '\alpha = 0.01';
res{ 7,1} = load('20150130-1531-SpinBoson-OrthPol-v37TCM33-alpha0.05delta0.1epsilon0dk20D5dopt5L200-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 7,2} = '\alpha = 0.05';
res{ 8,1} = load('20150130-1531-SpinBoson-OrthPol-v37TCM33-alpha0.1delta0.1epsilon0dk20D5dopt5L200-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 8,2} = '\alpha = 0.1';
res{ 9,1} = load('20150130-1531-SpinBoson-OrthPol-v37TCM33-alpha0.15delta0.1epsilon0dk20D5dopt5L200-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 9,2} = '\alpha = 0.15';
res{10,1} = load('20150130-1531-SpinBoson-OrthPol-v37TCM33-alpha0.2delta0.1epsilon0dk20D5dopt5L200-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{10,2} = '\alpha = 0.2';

%% TDVP SBM multi load files: Architecture sweep, OrthPol, artificial, L=50
res{ 1,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 1,2} = '\alpha = 0.01, Haswell';
res{ 2,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.05delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 2,2} = '\alpha = 0.05, Haswell';
res{ 3,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.1delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 3,2} = '\alpha = 0.1, Haswell';
res{ 4,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.15delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 4,2} = '\alpha = 0.15, Haswell';
res{ 5,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.2delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 5,2} = '\alpha = 0.2, Haswell';
res{ 6,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 6,2} = '\alpha = 0.01, Sandy';
res{ 7,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.05delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 7,2} = '\alpha = 0.05, Sandy';
res{ 8,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.1delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 8,2} = '\alpha = 0.1, Sandy';
res{ 9,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.15delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 9,2} = '\alpha = 0.15, Sandy';
res{10,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.2delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{10,2} = '\alpha = 0.2, Sandy';
res{11,1} = load('20150201-1747-SpinBoson-OrthPol-v37TCM14-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{11,2} = '\alpha = 0.01, Nehalem';
res{12,1} = load('20150201-1747-SpinBoson-OrthPol-v37TCM14-alpha0.05delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{12,2} = '\alpha = 0.05, Nehalem';
res{13,1} = load('20150201-1747-SpinBoson-OrthPol-v37TCM14-alpha0.1delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{13,2} = '\alpha = 0.1, Nehalem';
res{14,1} = load('20150201-1747-SpinBoson-OrthPol-v37TCM14-alpha0.15delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{14,2} = '\alpha = 0.15, Nehalem';
res{15,1} = load('20150201-1747-SpinBoson-OrthPol-v37TCM14-alpha0.2delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{15,2} = '\alpha = 0.2, Nehalem';
res{16,1} = load('20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpandBondExpand-small.mat');
res{16,2} = 'v19, Perfect';
%% TDVP SBM multi: Plot Visibility / Coherence
fignum = 2; figure(fignum); clf; hold all;
% pick = [1:length(res)];			% plot all
pick = [1,15:18];						% plot selective
plot(1:400,zeros(1,400),'black');
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
% ph{8}.LineStyle = ':';				%temp
set(gca,'ylim',[-1,1]);
xlabel('t');
ylabel('$<s_z>$');
legend([ph{:}],res{pick,2},'location','best');
formatPlot(fignum);

%% TDVP SBM multi: Plot VMPS GS <n>
fignum = 2; figure(fignum); clf; hold all;
pick = [1:length(res)];			% plot all
% pick = [8:11];
% pick = [23,35];
% pick = [6,12,21,22,23,24,25,26];						% plot selective
% pick = [8,15,16,11,12,13,14];
ph = cellfun(@(x) plot(real(x.results.nx)), res(pick,1), 'UniformOutput', false);
% ph{2}.LineStyle = ':';				%temp
% ph{4}.LineStyle = ':';				%temp
set(gca,'YScale','log');
xlabel('Site k')
ylabel('$<n_{k,VMPS}>$')
set(gca,'yscale','log')
legend([ph{:}],res{pick,2},'location','best');
formatPlot(fignum)

%% TDVP SBM multi: Plot GS Energy convergence
fignum = 3; figure(fignum); clf; hold all;
% pick = [1:length(res)];			% plot all
% pick = [8,15,16,11,12,13,14];						% plot selective
% pick = [8:11];
% pick = [23,35];
ph = cellfun(@(x) plot(cell2mat(x.results.EvaluesLog)-min(cell2mat(x.results.EvaluesLog))), res(pick,1), 'UniformOutput', false);
% disp(sprintf('%.15e\n',cell2mat(cellfun(@(x) x.results.E, res(pick,1), 'UniformOutput', false))))
set(gca,'YScale','log');
% try
% title(sprintf('$E_0 = %.10g, \\Lambda =  %.2g, z =  %.2g$',results.E, para.Lambda, para.z));catch end
xlabel('Site$\cdot$Loop');
ylabel('$E-E_0$');
legend([ph{:}],res{pick,2},'location','best');
formatPlot(fignum)
yLim = get(gca,'YLim');
for i = 1:para.loop
%     line([para.L*i para.L*i],yLim,'LineWidth',1,'Color','black');
end
if wantSave
    export_fig(sprintf('%s%s-MLSBM-Econvergence-Lambda%.2gz%.2gp16',saveto,para.filename(1:13),para.Lambda,para.z),'-transparent','-png','-painters')
end

%% TEST SBM_genpara:
para1.chain.mapping = 'OrthogonalPolynomials';	para1.chain.spectralDensity = 'Leggett_Hard'; para1.chain.discrMethod = 'Analytic';
para1.chain.discretization = 'None';			para1.chain.method = 'Analytic';
para1.Lambda = 2; para1.z = 1; para1.L = 100; para1.s = 1; para1.alpha = 0.1;
para1 = SBM_genpara(para1);
para2.chain.mapping = 'OrthogonalPolynomials';	para2.chain.spectralDensity = 'Leggett_Hard'; para2.chain.discrMethod = 'Numeric';
para2.chain.discretization = 'LogZ';			para2.chain.method = 'Stieltjes';
para2.Lambda = 2; para2.z = 1; para2.L = 100; para2.s = 1; para2.alpha = 0.1;
para2 = SBM_genpara(para2);
%%
% clear;
para1.chain.mapping = 'OrthogonalPolynomials';	para1.chain.spectralDensity = 'Leggett_Hard'; para1.chain.discrMethod = 'None';
para1.chain.discretization = 'LogZ';			para1.chain.method = 'Analytic';
para1.Lambda = 1.1; para1.z = 1; para1.L = 50; para1.s = 1; para1.alpha = 0.1;
t1 = cputime; para1 = SBM_genpara(para1); para1.time = cputime-t1;
para2.chain.mapping = 'OrthogonalPolynomials';	para2.chain.spectralDensity = 'Leggett_Hard'; para2.chain.discrMethod = 'Analytic';
para2.chain.discretization = 'LogZ';			para2.chain.method = 'Stieltjes';
para2.Lambda = 1.001; para2.z = 1; para2.L = 50; para2.s = 1; para2.alpha = 0.1;
t1 = cputime; para2 = SBM_genpara(para2); para2.time = cputime-t1;
%%
hold all;
% plot(para1.t);plot(para1.epsilon);
% plot(para2.t);plot(para2.epsilon);
plot(para1.t-para2.t);plot(para1.epsilon-para2.epsilon);
% set(gca,'yscale','log');
legend('para1.t','para1.\epsilon','para2.t',num2str(para2.Lambda));
% legend('t Lanzcos','\epsilon','t Stieltjes','\epsilon');
formatPlot(1)
