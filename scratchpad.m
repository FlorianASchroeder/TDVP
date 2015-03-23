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
i = 13;
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

%% Simulate shifting procedure
dim = 200; [bp,bm,n]=bosonop(dim,0,'n');
v1 = zeros(dim,1); v1(2:4) = 1/sqrt(3);
% 1. iteration
x = (bp+bm)./sqrt(2);
shift1 = v1'*x*v1;
[bp1n,bm1n,n1n]=bosonop(dim,shift1,'n');
bm1t = expm((bp-bm).*shift1./sqrt(2))'*bm*expm((bp-bm).*shift1./sqrt(2));
v1n = expm((bp-bm).*shift1./sqrt(2))'*v1;
% 2. iteration
x1n = (bp1n+bm1n)./sqrt(2);
shift2 = v1n'*x1n*v1n;
[bp2n,bm2n,n2n]=bosonop(dim,shift2,'n');
v2n = expm((bm1n-bp1n).*shift2./sqrt(2))*v1n;

plot([v1,v1n,v2n])
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
% set(gca,'yscale','log')
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

%% TDVP (1) SBM: Plot evolution of the spin
figure(1);clf;
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

%% TDVP (2) SBM: Plot <sz> / Visibility / Coherence
figure(2); hold all;
% plot(para.tdvp.t, tresults.spin.visibility);
plot(para.tdvp.t(1:length(tresults.spin.sz)), tresults.spin.sx);
plot(para.tdvp.t(1:length(tresults.spin.sz)), tresults.spin.sy);
plot(para.tdvp.t(1:length(tresults.spin.sz)), tresults.spin.sz);
set(gca,'ylim',[-1,1]);
% set(gca,'xscale','log');
xlabel('t');
ylabel('$\left<s_z\right>$');
legend('s_x','s_y','s_z');

%% TDVP (2) SBM: Plot Bloch length
figure(2); hold all;
plot(para.tdvp.t(1:length(tresults.spin.sz)), sqrt(tresults.spin.sz.^2+tresults.spin.sx.^2+tresults.spin.sy.^2));
plot(para.tdvp.t, tresults.spin.visibility);
set(gca,'ylim',[0,1]);
% set(gca,'xscale','log');
xlabel('t');
ylabel('$\sqrt{<s_x>^2+<s_y>^2+<s_z>^2}$');
legend('Bloch length','Visibility');

%% TDVP (3): Plot <n> environment chain
mode = 1;		% 0: lin, 1: log
figure(5); clf;
n = size(tresults.nx,1);
if isfield(tresults,'t')
	t=tresults.t;		% for the new convention when extracting in intervals >= rev42
else
	t=para.tdvp.t;		% for the old files
end
if mode
	surf(1:para.L,t(1:n),log10(real(tresults.nx)))
else
	surf(1:para.L,t(1:n),real(tresults.nx));
end
cb = colorbar;cb.Title.Interpreter = 'latex';
cb.Title.String = '$\left<n_k\right>$';
xlabel('Site $k$');
ylabel('Time $\omega_c t$');
zlabel('$<n_k>$');
% set(gca,'yscale','log');se
% set(gca,'View',[0 42]);
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
if mode
	set(gca,'zlim',[-4.5,0]);
	set(gca,'clim',[-4.5,0]);
	cb.Title.String = '$\log_{10}\left<n_k\right>$';
end

%% TDVP (3.1): Plot <n> star
mode = 1;		% 0: lin, 1: log
figure(6); clf;
if mode
% 	surf(tresults.star.omega,tresults.star.t,log10(tresults.star.n));
	surf(tresults.star.omega,tresults.star.t,log10(tresults.star.n)-log10(ones(size(tresults.star.n,1),1)*tresults.star.n(1,:)));

	zlabel('$\log_{10}\left<n_k\right>$');
else
	surf(tresults.star.omega,tresults.star.t,tresults.star.n);
	zlabel('$\left<n_k\right>$');
end
cb = colorbar;cb.Title.Interpreter = 'latex';
cb.Title.String = '$\left<n_k\right>$';
xlabel('Mode $\omega_k / \omega_c$');
ylabel('Time $\omega_c t$');
% set(gca,'yscale','log');		% log in sites
% set(gca,'View',[0 42]);
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
if mode
% 	set(gca,'zlim',[-4,1]);
% 	set(gca,'clim',[-4,1]);
	cb.Title.String = '$\log_{10}\left<n_k\right>$';
end

%% TDVP (3.2): Animate <n> propagation
figure(3); clf;
ax = axes('units','pixels');
pl = plot(1:para.L,log10(real(tresults.nx(1,:))));
set(gca,'ylimmode','manual');
set(gca,'ylim',[-5,0]);
set(gca,'xlimmode','manual','xlim',[1,para.L]);
xlabel('Site $k$');
ylabel('$\left<n_k\right>$');
% shading interp
sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',size(tresults.nx,1),'Value',1,...
        'Position', [400 20 120 20],...
        'Callback', @(source,callbackdata) set(pl,'ydata',log10(real(tresults.nx(round(source.Value),:)))));
%% TDVP (3.2): Animate <n> star formation
mode = 0;		% 0: lin, 1: log
figure(3); clf;
ax = axes('units','pixels');
if mode

	pl = plot(tresults.star.omega,log10(tresults.star.n(1,:)));
	ylabel('$log_{10}\left<n_k\right>$');
	set(gca,'ylim',[-5,log10(max(max(tresults.star.n)))]);
else
	pl = plot(tresults.star.omega,tresults.star.n(1,:));
	ylabel('$\left<n_k\right>$');
	set(gca,'ylim',[0,max(max(tresults.star.n))]);
end
set(gca,'ylimmode','manual');
set(gca,'xlimmode','manual','xlim',[0,tresults.star.omega(end)]);
xlabel('Site $k$');
% shading interp
if mode
	sld = uicontrol('Style', 'slider',...
			'Min',1,'Max',size(tresults.star.n,1),'Value',1,...
			'Position', [400 20 120 20],...
			'Callback', @(source,callbackdata) set(pl,'ydata',log10(tresults.star.n(round(source.Value),:))));
else
	sld = uicontrol('Style', 'slider',...
			'Min',1,'Max',size(tresults.star.n,1),'Value',1,...
			'Position', [400 20 120 20],...
			'Callback', @(source,callbackdata) set(pl,'ydata',tresults.star.n(round(source.Value),:)));
end
%% TDVP: Draw <n> propagation
figure(3); clf; hold on
n = size(tresults.nx,1); nCont = 30; picks = 8;
contourf(1:para.L,para.tdvp.t(1:n),real(tresults.nx),nCont)
%%
coneParam = zeros(2,picks);
hold on
for m=1:picks
	[coneParam(1,m),coneParam(2,m)] = ginput(1);
	plot(coneParam(1,m),coneParam(2,m),'ko','markerSize',10,'lineWidth',1,'markerFaceColor','w');
end
%%
% pf = polyfit(coneParam(1,:),coneParam(2,:),3);
% coneFunct = polyval(pf,para.tdvp.t(1:n));
pf = fnxtr(csaps(coneParam(2,:),coneParam(1,:)));
L = 800;
coneFunct = fnval(pf,1:L);
plot(coneFunct,1:L,'w','linewidth',4)
plot(coneFunct,1:L,'k','linewidth',2)

%% TDVP (4): Plot d_opt
figure(4); clf;
n = size(tresults.nx,1);
surf(1:para.L,para.tdvp.t(1:n),results.tdvp.d_opt)
xlabel('Site $k$');
ylabel('Time $t$');
zlabel('$d_{opt}$');
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight

%% TDVP (5): Plot D
figure(5); clf;
n = size(tresults.nx,1);
surf(1:para.L-1,para.tdvp.t(1:n),results.tdvp.D)
xlabel('Site $k$');
ylabel('Time $t$');
zlabel('$d_{opt}$');
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight

%% TDVP (6): Plot vNE of A / V
figure(7); clf;
results.tdvp.Amat_vNE = cell2mat(cellfun(@(x) sum(-x.^2.*log(x.^2)), results.tdvp.Amat_sv, 'UniformOutput',false));
results.tdvp.Vmat_vNE = cell2mat(cellfun(@(x) sum(-x.^2.*log(x.^2)), results.tdvp.Vmat_sv, 'UniformOutput',false));
subplot(1,2,1);
surf(1:size(results.tdvp.Amat_vNE,2),para.tdvp.t(1:size(results.tdvp.Amat_vNE,1)),results.tdvp.Amat_vNE);
xlabel('Bond $k$'); ylabel('Time $\omega_c t$'); zlabel('$S_{vNE}(A)$');
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
subplot(1,2,2);
surf(1:size(results.tdvp.Vmat_vNE,2),para.tdvp.t(1:size(results.tdvp.Vmat_vNE,1)),results.tdvp.Vmat_vNE);
xlabel('Bond $k$'); ylabel('Time $\omega_c t$'); zlabel('$S_{vNE}(V)$');
set(gca,'View',[0 90]);
shading interp
rotate3d on
axis tight
% formatPlot(7)

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

%% TDVP: Plot Min of Amat SV over time
figure(5);hold on;
minSV = cell2mat(cellfun(@min, results.tdvp.Amat_sv, 'UniformOutput',false));
surf(1:(para.L-1),para.tdvp.t(1:size(minSV,1)),log10(minSV))
% minSV = cell2mat(cellfun(@(x) log10(x(end-2:end)), results.Amat_sv(:,2:end),'UniformOutput',false));
% plot(minSV')
xlabel('Site $k$');
ylabel('Time $t$');
axis tight
shading interp
rotate3d on
set(gca,'zlim',[-10,0]);
%% TDVP: Plot Min of Vmat SV over time
figure(5);clf;
% minSV = cell2mat(cellfun(@(x) log10(min(x)), results.tdvp.Vmat_sv(:,2:end),'UniformOutput',false));
minSV = cell2mat(cellfun(@(x) log10(x(1:3)), results.Vmat_sv(:,2:end),'UniformOutput',false));
plot(minSV')
% surf(1:(para.L-1),para.tdvp.t(1:size(minSV,1)),minSV)
xlabel('Site $k$');
ylabel('Time $t$');
axis tight
shading interp
rotate3d on
% set(gca,'zlim',[-10,0]);

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
%%
pl = surf(plotMatTime{2,1});
title('Temporal change of OBB')
ylabel('$d_k$')
xlabel('Site $k$')
set(gca,'View',[9.5 40]);
formatPlot(6)
axis tight
rotate3d on
sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',size(plotMatTime,1),'Value',1,...
        'Position', [1100 20 120 20],...
        'Callback', @(source,callbackdata) set(pl,'zdata',plotMatTime{round(source.Value)}));

%% TDVP z-averaging in one file
% naming scheme to find files:
%   take series filename and replace z-value by *
clear;
folder = '20150225-1826-SpinBoson-LogZ-v40TCMde9-alpha1.5delta0.1Lambda1.1dk20D5dopt5L100-artificial';
figure(7);clf;
% folder = '20141117-0531-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L84';
% folder = '20141117-0406-SpinBoson-alpha0.2delta0.1epsilon0dk20D5dopt5L49';
filescheme = 'results-Till1000Step2v40-OBBExpand-noBondExpand-expvCustom800-1core-z*-small.mat';
files = ls(sprintf('%s/%s',folder,filescheme));
PlotData.spin.sz = [];PlotData.spin.sx = []; PlotData.spin.sy = []; PlotData.nx = [];
PlotData.z = [];
for k = 1:size(files,1)
    load([folder,'/',files(k,:)],'para','tresults');
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
	save(para.tdvp.filename, 'para','tresults');
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
set(gca,'xlim',[1,1e5]); set(gca,'xtick',[1,10,100,1e3,1e4,1e5]);
legend('Custom Krylov e^{At}v','Expokit')
xlabel('Matrix dimension n')
ylabel('Time/s')
formatPlot(9)

%% TDVP (12) expError analysis

figure(12); clf;
[n,m] = size(results.tdvp.expError);
surf(1:m,para.tdvp.t(2:end),log10(results.tdvp.expError))

xlabel('n-th Exp()');
ylabel('Time $t$');
zlabel('$log_{10}$(Exponential error)');
set(gca,'yscale','log');
% set(gca,'zscale','log');
set(gca,'View',[35 42]);
shading interp
rotate3d on
axis tight

%% TDVP (13) prepare artificial GroundState
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

%% TDVP SBM multi load files: OrthPol/LogZ v37: L=50, L=200; a = 0.01;  VMPS vs artificial		LabBook: 27/01/2015
% works only in cacheComputations!
defPlot(1,:) = {'OrthPol-TDVP-Benchmark-Sz-L50-v37',				[1,3:6]};
defPlot(2,:) = {'OrthPol-TDVP-Benchmark-Sz-L200-v37',				[1,7:9]};
defPlot(3,:) = {'OrthPol-TDVP-Benchmark-Sz-L50-artificial-v37',		[1,10:14]};
defPlot(4,:) = {'OrthPol-TDVP-Benchmark-Sz-L200-artificial-v37',	[1,15:18]};
defPlot(5,:) = {'OrthPol-TDVP-Benchmark-Sz-LogZ-Best-v37',			[1,2,11,16]};
% Legend font size: 16 or 18
cd('./../cacheComputations/');
res{1,1} = load('20141114-1902-SpinBoson-OrthPol-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4-OBBExpandBondExpand-small.mat','para','tresults');
res{1,2} = 'v19, Perfect';
res{2,1} = load('20150126-1800-SpinBoson-LogZ-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L49/results-Till325Step4v37-OBBExpand-noBondExpand-1core-zAvg.mat','para','tresults');
res{2,2} = '\Lambda=2, z_{Avg}=0.2, OBB no Bond';
res{3,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4v37-noOBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{3,2} = 'OrthPol L=50, no Expand';
res{4,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4v37-OBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{4,2} = 'OrthPol L=50, OBB no Bond';
res{5,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4v37-noOBBExpand-BondExpand15-1core-small.mat','para','tresults');
res{5,2} = 'OrthPol L=50, Bond no OBB';
res{6,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom800-1core-small.mat','para','tresults');
res{6,2} = 'OrthPol L=50, All Expand';
res{7,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v37-noOBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{7,2} = 'OrthPol L=200, no Expand';
res{8,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v37-OBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{8,2} = 'OrthPol L=200, OBB no Bond';
res{9,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom800-1core-small.mat','para','tresults');
res{9,2} = 'OrthPol L=200, All Expand';
res{10,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-noOBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{10,2} = 'OrthPol L=50, Art, no Expand';
res{11,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-OBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{11,2} = 'OrthPol L=50, Art, OBB no Bond';
res{12,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-noOBBExpand-BondExpand15-1core-small.mat','para','tresults');
res{12,2} = 'OrthPol L=50, Art, Bond no OBB';
res{13,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom800-1core-small.mat','para','tresults');
res{13,2} = 'OrthPol L=50, Art, All, expvCustom800';
res{14,1} = load('20150126-1718-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom1-1core-small.mat','para','tresults');
res{14,2} = 'OrthPol L=50, Art, All, expvCustom1';
res{15,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial\results-Till325Step4v37-noOBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{15,2} = 'OrthPol L=200,Art, no Expand';
res{16,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial\results-Till325Step4v37-OBBExpand-noBondExpand-1core-small.mat','para','tresults');
res{16,2} = 'OrthPol L=200,Art, OBB no Bond';
res{17,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial\results-Till325Step4v37-noOBBExpand-BondExpand15-1core-small.mat','para','tresults');
res{17,2} = 'OrthPol L=200,Art, Bond no OBB';
res{18,1} = load('20150126-1719-SpinBoson-OrthPol-v37TCM68-alpha0.01delta0.1epsilon0dk20D5dopt5L200-artificial\results-Till325Step4v37-OBBExpand-BondExpand15-expvCustom800-1core-small.mat','para','tresults');
res{18,2} = 'OrthPol L=200,Art, All';
res{19,1} = load('20150126-2247-SpinBoson-LogZ-v37TCM66-alpha0.2delta0.1epsilon0dk20D5dopt5L49/results-Till325Step4v37-OBBExpand-noBondExpand-1core-zAvg.mat','para','tresults');
res{19,2} = '\Lambda=2, z_{Avg}=0.2, OBB no Bond';
res{20,1} = load('20150126-2247-SpinBoson-LogZ-v37TCM66-alpha0.2delta0.1epsilon0dk20D5dopt5L49-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-1core-zAvg.mat','para','tresults');
res{20,2} = '\Lambda=2, z_{Avg}=0.2, OBB no Bond, Art';

for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight
	% ph{8}.LineStyle = ':';				%temp
	set(gca,'ylim',[-1,1]);
	xlabel('t');
	ylabel('$<s_z>$');
	leg = legend([ph{:}],res{pick,2},'location','best');
	leg.FontSize = 18;
	formatPlot(fignum);
	title(defPlot{fignum,1},'fontsize',15);		% comment me
end
cd('./../TDVP/');
%% TDVP SBM multi load files: Architecture sweep, OrthPol, artificial, L=50						LabBook: 01/02/2015
col = get(groot,'defaultAxesColorOrder');
defPlot(1,:) = {'Orth2010-OrthPol-TDVP-OBBnoBondExpand-L50-architecture-artificial-v37',	[1:15]};

res{ 1,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 1,2} = '\alpha = 0.01, Haswell'; res{1,3} = 'Haswell';
res{ 2,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.05delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 2,2} = '\alpha = 0.05, Haswell';
res{ 3,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.1delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 3,2} = '\alpha = 0.1, Haswell';
res{ 4,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.15delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 4,2} = '\alpha = 0.15, Haswell';
res{ 5,1} = load('20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha0.2delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 5,2} = '\alpha = 0.2, Haswell';
res{ 6,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 6,2} = '\alpha = 0.01, Sandy'; res{6,3} = 'Sandy';
res{ 7,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.05delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 7,2} = '\alpha = 0.05, Sandy';
res{ 8,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.1delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 8,2} = '\alpha = 0.1, Sandy';
res{ 9,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.15delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{ 9,2} = '\alpha = 0.15, Sandy';
res{10,1} = load('20150201-1745-SpinBoson-OrthPol-v37TCM40-alpha0.2delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{10,2} = '\alpha = 0.2, Sandy';
res{11,1} = load('20150201-1747-SpinBoson-OrthPol-v37TCM14-alpha0.01delta0.1epsilon0dk20D5dopt5L50-artificial/results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{11,2} = '\alpha = 0.01, Nehalem'; res{11,3} = 'Nehalem';
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
for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');

	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight
	cellfun(@(x) set(x,'Color',col(1,:)), ph(1:5));
	cellfun(@(x) set(x,'Color',col(2,:)), ph(6:10));
	cellfun(@(x) set(x,'Color',col(3,:)), ph(11:15));
	set(gca,'ylim',[-1,1]);
	xlabel('t');
	ylabel('$<s_z>$');
	leg = legend([ph{[1,6,11]}],res{[1,6,11],3},'location','best');
	leg.FontSize = 18;
	formatPlot(fignum);
	title(defPlot{fignum,1},'fontsize',15);
end
%% TDVP SBM Ohmic  s1  0.5>a  Orth2010, OrthPol, art, L=50, L=200; weak -0.75 coupling		LabBook: 30/01/2015, 27/02/2015, POSTER April 2015
cd('./../cacheComputations/');
defPlot(1,:) = {'Orth2010-OrthPol-TDVP-OBBnoBondExpand-L50-artificial-v37',					[1:5],			{'ylim',[-1,1]}};
defPlot(2,:) = {'Orth2010-OrthPol-TDVP-OBBnoBondExpand-L200-artificial-v37',				[6:10],			{'ylim',[-1,1]}};
defPlot(3,:) = {'Orth2010-OrthPol-TDVP-s1-med-AllExpand-L200-artificial-v41',				[11:12, 20:21],	{'ylim',[-0.2,1],'xlim',[0,320]}};		% POSTER, Incomplete, PC67
% defPlot(4,:) = {'Orth2010-13a-OrthPol-TDVP-OBBExpand-L200-DeltaT1-artificial-v40',			[15:19],		{'ylim',[-1,1]}};		% This is the PERFECT one!
defPlot(4,:) = {'Orth2010-13a-OrthPol-TDVP-OBBExpand-L200-DeltaT1-artificial-v40',			[15:19],		{'ylim',[-1,1],'xlim',[0,320]}};		% POSTER, This is the PERFECT one!

% 01-05
foldPattern = '20150130-1334-SpinBoson-OrthPol-v37TCM33-alpha*delta0.1epsilon0dk20D5dopt5L50-artificial';
filePattern = 'results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = cell(size(folds,1),4);
i = 0;
for file = {folds.name}
	i = i+1;
	file = file{1};
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res = sortrows(res,3); res{1,4} = 'OBBExpand, L=50, art.';

% 06-10
foldPattern = '20150130-1531-SpinBoson-OrthPol-v37TCM33-alpha*delta0.1epsilon0dk20D5dopt5L200-artificial';
filePattern = 'results-Till325Step4v37-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, L=200, art.';
cd('./../TDVP/');

% 11-14		Medium			Finer REDO!! look at PC67
foldPattern = '20150218-1340-SpinBoson-OrthPol-v39TCM1-alpha*delta0.1epsilon0dk20D5dopt5L200-artificial';
filePattern = 'results-Till325Step0.5v41-OBBExpand-BondExpand7-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, L=200, art.';

% 15-19
foldPattern = '20150227-231*-SpinBoson-OrthPol-v40TCMde9-alpha*delta0.1Lambda2dk20D5dopt5L200-artificial';
filePattern = 'results-Till1000Step1v40-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, L=200, art., T=1000step1';

% 20-21
foldPattern = '20150305-0328-SpinBoson-OrthPol-v41TCMde9-alpha*delta0.1epsilon0dk20D5dopt5L400-artificial';
filePattern = 'results-Till2000Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, L=400, t=0.5';


for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],res{pick,2},'location','best');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum);
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=1$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
% 	title(defPlot{fignum,1},'fontsize',15);
end

%% TDVP SBM Ohmic  s1  0.5<a  Orth2010, OrthPol, artificial, L=200; strong coupling
clear;
defPlot(1,:) = {'Orth2010-13b-OrthPol-TDVP-noExpand-L200-artificial-v40',						[1:5],	{'ylim',[-1,1]}};
defPlot(2,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L200-artificial-v40-expTime',				[6:10],	{'ylim',[-1,1]}};
defPlot(3,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L200-DeltaT2-artificial-v40',				[11:16],{'ylim',[-1,1]}};		% Incomplete?
defPlot(4,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L200-0.5-0.75-DeltaT0.5-artificial-v40',   [17:18,21:22],{'ylim',[0,1], 'xscale','log'}};		% Incomplete pc67
defPlot(5,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L200-DeltaT0.5-2-artificial-v40',			[17:18,14:16],{'ylim',[0,1], 'xscale','log','xlim',[1,1000]}};		% Incomplete pc67
defPlot(6,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L500-DeltaT0.01-v42',						[23:28],{'ylim',[-1,1], 'xscale','log','xlim',[1,320]}};		% Incomplete pc67

% 1-5
foldPattern = '20150225-1301*OrthPol*'; filePattern = '*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = cell(size(folds,1),4);
i = 0;
for file = {folds.name}
	i = i+1;
	file = file{1};
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res = sortrows(res,3); res{1,4} = 'noExpand';

i = i+1;
res{i,1} = load('20150219-0043-SpinBoson-OrthPol-v39TCM67-alpha0.5delta0.1epsilon0dk20D5dopt5L200-artificial/results-Till100000ExpStep0.1v39-OBBExpand-noBondExpand-expvCustom800-1core-small.mat','para','tresults');
res{i,2} = '\alpha = 0.5'; res{i,3} = 0.5;

% 7-10
foldPattern = '20150218-1609*OrthPol*'; filePattern = '*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-ExpTime';

% 11-16
foldPattern = '20150227-18*OrthPol*'; filePattern = '*Step2*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);

% 17-18
foldPattern = '20150305-0328-SpinBoson-OrthPol-v41TCMde9-alpha*delta0.1epsilon0dk20D5dopt5L400-artificial';
filePattern = 'results-Till2000Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-0.5new';

% 19-20
foldPattern = '20150227-18*OrthPol*'; filePattern = '*Step1*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-1';

% 21-22
foldPattern = '20150227-18*OrthPol*'; filePattern = '*Step0.5*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-0.5';

for fignum = 3:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight
	% ph{8}.LineStyle = ':';				%temp
	set(gca,defPlot{fignum,3}{:});
	xlabel('t');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],res{pick,2},'location','best');
	leg.FontSize = 18;
	formatPlot(fignum);
	set(gca,'color','none');
% 	title(defPlot{fignum,1},'fontsize',15);
end

%% TDVP SBM Ohmic  s1  0.4<a  Orth2010, OrthPol, coupled, L=500, NEW							LabBook: 08/05/15, Paper
clear;
defPlot(1,:) = {'Orth2010-13b-OrthPol-TDVP-OBBExpand-L500-DeltaT0.01-coupled-v42',						[1:6],{'ylim',[-0.1,1], 'xscale','log','xlim',[1,400]}};		% Incomplete pc67

res = cell(0,4);
i = 0;

foldPattern = '20150409-1304-SpinBoson-OrthPol-v42TCMde9-alpha*delta0.1epsilon0dk30D5dopt5L500';
filePattern = 'results-Till2000Step0.01v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand-dt0.01';

for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.tresults.t(sum(x.tresults.t~=0)), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.tresults.t(1:sum(x.tresults.t~=0)), x.tresults.spin.sz(1:sum(x.tresults.t~=0))), res(pick,1), 'UniformOutput', false);
	axis tight, ax = gca;
	% ph{8}.LineStyle = ':';				%temp
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum);
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=1$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
% 	title(defPlot{fignum,1},'fontsize',15);
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM subOhmic s0.5 SBM: Orth2010, OrthPol, L=200:										LabBook: 06/03/2015, Poster!
clear
defPlot( 1,:) = {'Orth2010-14ab-OrthPol-TDVP-noExpand-L200-DeltaT4-artificial-v40',		[1:9],	 {'ylim',[-1,1]}};
defPlot( 2,:) = {'Orth2010-14ab-OrthPol-TDVP-OBBExpand-L200-DeltaT4-artificial-v40',	[10:18], {'ylim',[-1,1]}};	% not good enough!
defPlot( 3,:) = {'Orth2010-14ab-OrthPol-TDVP-OBBExpand-L200-DeltaT0.5-artificial-v40',	[19:27], {'ylim',[-1,1]}};
defPlot( 4,:) = {'Orth2010-14ab-OrthPol-TDVP-OBBExpand-L100-DeltaT0.5-coupled-v41',		[28:36], {'ylim',[-1,1]}};
defPlot( 5,:) = {'Orth2010-14ab-OrthPol-TDVP-OBBExpand-L100-DeltaT0.5-dk50-coupled-v41',[37:45], {'ylim',[-1,1]}};
defPlot( 6,:) = {'Orth2010-14ab-OrthPol-TDVP-AllExpand-L100-DeltaT0.5-dk50-coupled-v41',[46:54], {'ylim',[-1,1]}};
defPlot( 7,:) = {'Orth2010-14a-OrthPol-TDVP-OBBExpand-L200-DeltaT0.5-artificial-v40',	[19:22], {'ylim',[-1,1]}};
defPlot( 8,:) = {'Orth2010-14b-OrthPol-TDVP-OBBExpand-L200-DeltaT0.5-artificial-v40',	[23:27], {'ylim',[-0.5,1]}};
defPlot( 9,:) = {'Orth2010-14a-OrthPol-TDVP-OBBExpand-L100-DeltaT0.5-coupled-v41',		[28:31], {'ylim',[-1,1],'xlim',[0,320]}};		% Poster
defPlot(10,:) = {'Orth2010-14b-OrthPol-TDVP-OBBExpand-L100-DeltaT0.5-coupled-v41',		[32:36], {'ylim',[0,1],'xlim',[0,320]}};		% Poster

% 01-09
foldPattern = '20150225-171*-SpinBoson-OrthPol-v40TCMde9-s0.5-alpha*delta0.1epsilon0dk20D5dopt5L200-artificial';
filePattern = '*Till324ExpStep4v40-noOBBExpand-noBondExpand*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = cell(size(folds,1),4);
i = 0;
for file = {folds.name}
	i = i+1;
	file = file{1};
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res = sortrows(res,3); res{1,4} = 'noExpand';

% 10-18
foldPattern = '20150225-171*-SpinBoson-OrthPol-v40TCMde9-s0.5-alpha*delta0.1epsilon0dk20D5dopt5L200-artificial';
filePattern = '*Till324ExpStep4v40-OBBExpand-noBondExpand*-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand';

% 19-27
foldPattern = '20150225-171*-SpinBoson-OrthPol-v40TCMde9-s0.5-alpha*delta0.1epsilon0dk20D5dopt5L200-artificial';
filePattern = '*Till500Step0.5v*-OBBExpand-noBondExpand*small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, \Deltat=0.5';

% 28-36
foldPattern = '20150311-172*-SpinBoson-OrthPol-v41TCMde9-s0.5-alpha*delta0.1epsilon0dk20D5dopt5L100';
filePattern = 'results-Till324Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, \Deltat=0.5';

% 37-45
foldPattern = '20150312-1059-SpinBoson-OrthPol-v41TCMde9-s0.5-alpha*delta0.1epsilon0dk50D5dopt5L100';
filePattern = 'results-Till324Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'OBBExpand, \Deltat=0.5, dk50';

% 46-54
foldPattern = '20150312-1059-SpinBoson-OrthPol-v41TCMde9-s0.5-alpha*delta0.1epsilon0dk50D5dopt5L100';
filePattern = 'results-Till324Step0.5v41-OBBExpand-BondExpand10-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,4} = 'AllExpand, \Deltat=0.5';

for fignum = 5:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum);
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
% 	title(defPlot{fignum,1},'fontsize',15);
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM subOhmic s0.25, s0.75 SBM: Kast2013, OrthPol, L=50:									LabBook, Poster!
% use: export_fig(defPlot{get(gcf,'number'),1},'-transparent','-png')
clear
defPlot(1,:) = {'Kast2013-2-OrthPol-TDVP-OBBExpand-L50-DeltaT0.5-artificial-v41',		[1:8],	 {'ylim',[-0.5,1]}};
defPlot(2,:) = {'Kast2013-2-OrthPol-TDVP-AllExpand-L50-DeltaT0.5-artificial-v41',		[9:16],  {'ylim',[-0.5,1]}};
defPlot(3,:) = {'Kast2013-2-OrthPol-TDVP-OBBExpand-L50-DeltaT0.5-coupled-v41',			[17:24], {'ylim',[-0.5,1]}};
defPlot(4,:) = {'Kast2013-3-OrthPol-TDVP-OBBExpand-L50-DeltaT0.5-artificial-v41',		[25:31], {'ylim',[-0.5,1]}};
defPlot(5,:) = {'Kast2013-3-OrthPol-TDVP-OBBExpand-L50-DeltaT0.5-coupled-v41',			[32:38], {'ylim',[-0.5,1]}};
defPlot(6,:) = {'Kast2013-OrthPol-TDVP-OBBExpand-L50-200-coupled-v41',					[39:41], {'ylim',[-0.1,1]}};

foldPattern = '20150307-0341-SpinBoson-OrthPol-v41TCMde9-s0.75-alpha*delta0.1epsilon0dk20D5dopt5L50-artificial';
filePattern = 'results-Till100Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = cell(size(folds,1),3);
i = 0;
for file = {folds.name}
	i = i+1;
	file = file{1};
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res = sortrows(res,3); res{1,3} = 'artificial';

foldPattern = '20150307-0341-SpinBoson-OrthPol-v41TCMde9-s0.75-alpha*delta0.1epsilon0dk20D5dopt5L50-artificial';
filePattern = 'results-Till100Step0.5v41-OBBExpand-BondExpand15-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),3)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,3} = 'AllExpand';

foldPattern = '20150307-1406-SpinBoson-OrthPol-v41TCMde9-s0.75-alpha*delta0.1epsilon0dk20D5dopt5L50';
filePattern = 'results-Till100Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),3)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,3} = 'coupled';

foldPattern = '20150307-1324-SpinBoson-OrthPol-v41TCMde9-s0.25-alpha*delta0.1epsilon0dk20D5dopt5L50-artificial';
filePattern = 'results-Till100Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),3)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,3} = 's0.25 artificial';

foldPattern = '20150307-1350-SpinBoson-OrthPol-v41TCMde9-s0.25-alpha*delta0.1epsilon0dk20D5dopt5L50';
filePattern = 'results-Till100Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),3)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,3} = 's0.25 coupled';

foldPattern = '20150307-154*-SpinBoson-OrthPol-v41TCMde9-s0.75-alpha0.14delta0.1epsilon0dk20D5dopt5L*';
filePattern = 'results-Till100Step0.5v41-OBBExpand-noBondExpand-expvCustom800-1core.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),3)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'opt5L')+5:strfind(file,'\')-1));
	res{i,2} = sprintf('L = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res{offset+1,3} = 's0.25 Longer';


for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight
	set(gca,defPlot{fignum,3}{:});
	% set(gca,'xlim',[1,1e8]);set(gca,'xscale','log');
	xlabel('t');
	ylabel('$<s_z>$');
	leg = legend([ph{:}],res{pick,2},'location','best');
	leg.FontSize = 18;
	set(gca,'color','none');
	formatPlot(fignum);
% 	title(defPlot{fignum,1},'fontsize',15);
	if fignum == 5
		ax = axes('Position',[0.17,0.225,0.4,0.3275]);hold on;
		ax.ColorOrderIndex=4; ax.FontSize = 14;
		ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz,...
			'LineWidth',2), res(pick(4:end),1), 'UniformOutput', false);
		box(gca,'on');
		set(gca,'ylim',[0.93,1],'xlim',[0,30]);

	end
end

% figure(fignum+1);clf; hold all;
% ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM compare ML-MCTDH
clear
defPlot(1,:) = {'Wang2010-4-OrthPol-TDVP-OBBExpand-L50-DeltaT0.5-artificial-v41',		1:3,			{'ylim',[-1,1], 'xlim',[0,120]}};
defPlot(2,:) = {'Wang2010-7b-OrthPol-TDVP-AllExpand-L200-Trial-artificial-v42',			[2,4,6,8,11],		{'ylim',[0.6,1],  'xlim',[0,40]}};
defPlot(3,:) = {'Wang2010-7c-OrthPol-TDVP-AllExpand-L200-Trial-artificial-v42',			[3,5,7,9],		{'ylim',[0.6,1], 'xlim',[0,20]}};
defPlot(4,:) = {'Wang2010-7a-OrthPol-TDVP-AllExpand-L600-Trial-artificial-v42',			[10],			{'ylim',[0.4,1], 'xlim',[0,20]}};

res = {};
res{1,1} = load('20150316-1346-SpinBoson-OrthPol-v41TCM51-s0.5-alpha0.005delta0.2epsilon0dk20D5dopt5L300-artificial/results-Till1200Step1v41-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{1,2} = '\alpha = 0.005, Fig4';				res{1,3} = '\Deltat=1,OBB';
res{2,1} = load('20150316-1347-SpinBoson-OrthPol-v41TCM51-s0.5-alpha0.5delta0.2epsilon0dk20D5dopt5L100-artificial/results-Till400Step0.5v41-OBBExpand-BondExpand7-expvCustom800-1core-small.mat');
res{2,2} = '\alpha = 0.5, Fig7b, \Deltat=0.5';	res{2,3} = '\Deltat=0.5,All7';
res{3,1} = load('20150316-1348-SpinBoson-OrthPol-v41TCM51-s0.5-alpha1delta0.2epsilon0dk20D5dopt5L100-artificial/results-Till200Step0.5v41-OBBExpand-BondExpand7-expvCustom800-1core-small.mat');
res{3,2} = '\alpha = 1, Fig7c, \Deltat=0.5';	res{3,3} = '\Deltat=0.5,All7';
res{4,1} = load('20150323-1555-SpinBoson-OrthPol-v42TCM66-s0.5-alpha0.5delta0.2epsilon0dk20D5dopt5L200-artificial/results-Till20Step0.1v42-OBBExpand-BondExpand20-expvCustom800-1core-small.mat');
res{4,2} = '\alpha = 0.5, Fig7b, \Deltat=0.1';	res{4,3} = '\Deltat=0.1,All20';
res{5,1} = load('20150323-1555-SpinBoson-OrthPol-v42TCM66-s0.5-alpha1delta0.2epsilon0dk20D5dopt5L200-artificial/results-Till20Step0.1v42-OBBExpand-BondExpand20-expvCustom800-1core-small.mat');
res{5,2} = '\alpha = 1, Fig7c, \Deltat=0.1';	res{5,3} = '\Deltat=0.1,All20';
res{6,1} = load('20150323-1555-SpinBoson-OrthPol-v42TCM66-s0.5-alpha0.5delta0.2epsilon0dk20D5dopt5L200-artificial/results-Till50Step0.02v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{6,2} = '\alpha = 0.5, Fig7b, \Deltat=0.02'; res{6,3} = '\Deltat=0.02,OBB';
res{7,1} = load('20150323-1555-SpinBoson-OrthPol-v42TCM66-s0.5-alpha1delta0.2epsilon0dk20D5dopt5L200-artificial/results-Till50Step0.02v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{7,2} = '\alpha = 1, Fig7c, \Deltat=0.02';	res{7,3} = '\Deltat=0.02,OBB';
res{8,1} = load('20150323-1555-SpinBoson-OrthPol-v42TCM66-s0.5-alpha0.5delta0.2epsilon0dk20D5dopt5L200-artificial/results-Till50Step0.02v42-OBBExpand-BondExpand20-expvCustom800-1core-small.mat');
res{8,2} = '\alpha = 0.5, Fig7b, \Deltat=0.02';	res{8,3} = '\Deltat=0.02,All20';
res{9,1} = load('20150323-1555-SpinBoson-OrthPol-v42TCM66-s0.5-alpha1delta0.2epsilon0dk20D5dopt5L200-artificial/results-Till50Step0.02v42-OBBExpand-BondExpand20-expvCustom800-1core-small.mat');
res{9,2} = '\alpha = 1, Fig7c, \Deltat=0.02';	res{9,3} = '\Deltat=0.02,All20';
res{10,1} = load('20150325-1532-SpinBoson-OrthPol-v42TCM59-s0.5-alpha0.3delta0.2epsilon0dk20D5dopt5L600-artificial/results-Till50Step0.02v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{10,2} = '\alpha = 0.3, Fig7a, \Deltat=0.02';	res{10,3} = '\Deltat=0.02,OBB';
res{11,1} = load('20150325-1532-SpinBoson-OrthPol-v42TCM59-s0.5-alpha0.5delta0.2epsilon0dk20D5dopt5L600-artificial/results-Till50Step0.02v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat');
res{11,2} = '\alpha = 0.5, Fig7b, \Deltat=0.02';	res{11,3} = '\Deltat=0.02,OBB';

for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)).*(-x.para.hx)./2, x.tresults.spin.sz,'*'), res(pick,1), 'UniformOutput', false);
	axis tight
	set(gca,defPlot{fignum,3}{:});
	xlabel('$t\Delta$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],res{pick,3},'location','best');
	leg.FontSize = 18;
	set(gca,'color','none');
	formatPlot(fignum);
% 	title(defPlot{fignum,1},'fontsize',15);
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM multi: Ohmic LogZ																	LabBook: 04/02/2015
clear
defPlot(1,:) = {'Orth2010-13b-LogZ1.1-TDVP-OBBExpand-L100-DeltaT2-artificial-v40',		[1:6]};

foldPattern = '20150225*LogZ*'; filePattern = '*Till1000Step2v40-OBBExpand-noBondExpand-expvCustom800-1core*-zAvg.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = cell(size(folds,1),3);
i = 0;
for file = {folds.name}
	i = i+1;
	file = file{1};
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('\\alpha = %g', res{i,3});
end
res = sortrows(res,3);

for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight
	set(gca,'ylim',[-1,1]);
% 	set(gca,'xlim',[1,1e3]);set(gca,'xscale','log');
	xlabel('t');
	ylabel('$<s_z>$');
	leg = legend([ph{:}],res{pick,2},'location','best');
	leg.FontSize = 18;
	formatPlot(fignum);
	title(defPlot{fignum,1},'fontsize',15);
end

%% TDVP SBM multi files: Star Bath occupation
clear
defPlot(1,:) = {'Orth2010-13a-s1-OrthPol-TDVP-OBB10-L200-325T1-art-v42',		[1:5],	{'ylim',[-1,1]}};
defPlot(2,:) = {'Orth2010-14a-s0.5-OrthPol-TDVP-OBB10-L100-325T1-art-v42',		[6:9],	{'ylim',[-1,1]}};
defPlot(3,:) = {'Orth2010-14b-s0.5-OrthPol-TDVP-OBB10-L100-325T1-art-v42',		[10:14],{'ylim',[-0.5,1]}};
defPlot(4,:) = {'Orth2010-13a-s1-OrthPol-TDVP-OBB10-L200-325T1-coup-v42',		[15:19],{'ylim',[-1,1]}};
defPlot(5,:) = {'Orth2010-14a-s0.5-OrthPol-TDVP-OBB10-L100-325T1-coup-v42',		[20:23],{'ylim',[-1,1]}};
defPlot(6,:) = {'Orth2010-14b-s0.5-OrthPol-TDVP-OBB10-L100-325T1-coup-v42',		[24:28],{'ylim',[0,1]}};


foldPattern = '20150327-1327-SpinBoson-OrthPol-v42TCMde9-alpha*delta0.1epsilon0dk20D5dopt5L200-artificial';
filePattern = 'results-Till325Step1v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = cell(size(folds,1),4);
i = 0;
for file = {folds.name}
	i = i+1;
	file = file{1};
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res = sortrows(res,3); res(1:i,4) = {'$s=1$, vac'};

foldPattern = '20150327-1330-SpinBoson-OrthPol-v42TCMde9-s0.5-alpha*delta0.1epsilon0dk20D5dopt5L100-artificial';
filePattern = 'results-Till325Step1v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,4) = {'$s=0.5$, vac'};

foldPattern = '20150327-1430-SpinBoson-OrthPol-v42TCMde10-alpha*delta0.1epsilon0dk20D5dopt5L100';
filePattern = 'results-Till325Step1v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,4) = {'$s=1$, coup'};

foldPattern = '20150327-1434-SpinBoson-OrthPol-v42TCMde10-s0.5-alpha*delta0.1epsilon0dk20D5dopt5L100';
filePattern = 'results-Till325Step1v42-OBBExpand-noBondExpand-expvCustom800-1core-small.mat';
folds = rdir([foldPattern,'\',filePattern]);
res = [res;cell(size(folds,1),4)]; offset = i;
for file = {folds.name}
	file = file{1};
	i = i+1;
	res{i,1} = load(file,'para','tresults');			% comment first!
	res{i,3} = str2double(file(strfind(file,'alpha')+5:strfind(file,'delta')-1));
	res{i,2} = sprintf('$\\alpha$ = %g', res{i,3});
end
res((offset+1):end,:) = sortrows(res((offset+1):end,:),3);
res(offset+1:i,4) = {'$s=0.5$, coup'};

for fignum = 1:size(defPlot,1)
	f = figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight; ax = gca;
	set(ax,defPlot{fignum,3}{:});
	xlabel('$\omega_ct$');
	ylabel('$\left<\sigma_z\right>$');
	leg = legend([ph{:}],cellfun(@(x) sprintf('%g\n',x),res(pick,3),'UniformOutput',false),'location','best');
	legend boxoff
	fs = 22;
	leg.FontSize = fs;
	formatPlot(fignum);
	t1 = text(leg.Position(1)+ax.Position(1),leg.Position(2)+leg.Position(4)/2,'$\alpha$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	t2 = text(leg.Position(1),leg.Position(2)+leg.Position(4)/2,'$s=0.5$', 'FontSize',fs,'Units','norm','VerticalAlignment','bottom');
	set(gca,'color','none');
% 	title(defPlot{fignum,1},'fontsize',15);
end

figure(fignum+1);clf; hold all;
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.para.tdvp.calcTime)), x.para.tdvp.calcTime), res(pick,1), 'UniformOutput', false);

%% TDVP SBM multi (1): Plot Visibility / Coherence
fignum = 3; figure(fignum); clf; hold all;
pick = [1:length(res)];			% plot all
% pick = [1:5];						% plot selective
xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
plot([1,xmax],[0,0],'black');
ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
% ph{8}.LineStyle = ':';				%temp
set(gca,'ylim',[-1,1]);
% set(gca,'xlim',[0,325]);
% set(gca,'xlim',[1,xmax]);set(gca,'xscale','log');
xlabel('t');
ylabel('$\left<\sigma_z\right>$');
legend([ph{:}],res{pick,2},'location','best');
formatPlot(fignum);

%% TDVP SBM multi (2): Plot Visibility / Coherence		Predefined
for fignum = 1:size(defPlot,1)
	figure(fignum); clf; hold all;
	pick = defPlot{fignum,2};			% plot all
	xmax = max(cellfun(@(x) x.para.tdvp.t(end), res(pick,1)));
	plot([1,xmax],[0,0],'black');
	ph = cellfun(@(x) plot(x.para.tdvp.t(1:length(x.tresults.spin.sz)), x.tresults.spin.sz), res(pick,1), 'UniformOutput', false);
	axis tight
% 	cellfun(@(x) set(x,'Color',col(1,:)), ph(1:5));
% 	cellfun(@(x) set(x,'Color',col(2,:)), ph(6:10));
% 	cellfun(@(x) set(x,'Color',col(3,:)), ph(11:15));
% 	ph{mod(1:15,3)}. = ':';				%temp
	set(gca,'ylim',[-1,1]);
	% set(gca,'xlim',[1,1e8]);set(gca,'xscale','log');
	xlabel('t');
	ylabel('$<s_z>$');
	leg = legend([ph{:}],res{pick,2},'location','best');
	leg.FontSize = 18;
	formatPlot(fignum);
	title(defPlot{fignum,1},'fontsize',15);
end

%% TDVP SBM multi (3): Plot VMPS GS <n>
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

%% TDVP SBM multi (4): Plot GS Energy convergence
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

%% TDVP SBM multi (5): plot single file STAR
mode = 1;		% 0: lin, 1: log
fignum = 7;
figure(fignum); clf; x = res(1,:);
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1),pos(2),840,pos(4)]);pos = get(gcf,'Position');
pl = surf(log10(x{1}.tresults.star.n));
pl.UserData=x;
cb = colorbar;cb.Title.Interpreter = 'latex';
if mode
	zlabel('$log_{10}\left<n_k\right>$');
	pl.ZDataSource = 'log10(pl.UserData{1}.tresults.star.n);';
% 	pl.ZDataSource = 'log10(pl.UserData{1}.tresults.star.n)-log10(ones(size(pl.UserData{1}.tresults.star.n,1),1)*pl.UserData{1}.tresults.star.n(1,:));';
% 	set(cb,'defaulttextinterpreter','tex')
	cb.Title.String = 'log$_{10}\left<n_k\right>$';
else
% 	pl = surf(x{1}.tresults.star.omega,x{1}.tresults.star.t,x{1}.tresults.star.n);
	zlabel('$\left<n_k\right>$');
	pl.UserData=x;
	pl.ZDataSource = 'pl.UserData{1}.tresults.star.n;';
	cb.Title.String = '$\left<n_k\right>$';
end
pl.XDataSource = 'pl.UserData{1}.tresults.star.omega';
pl.YDataSource = 'pl.UserData{1}.tresults.star.t';refreshdata;
xlabel('Mode $\omega_k / \omega_c$');
ylabel('Time $\omega_c t$');
% set(gca,'zscale','log');
set(gca,'View',[0 90],'TickDir','out','FontSize',14);
shading interp;rotate3d on;axis tight;
if mode
% 	set(gca,'zlim',[0,10]);
% 	set(gca,'clim',[0,10]);
end
set(gca,'xlimmode','manual','xlim',[0,x{1}.tresults.star.omega(end)]);
com = ', ';
% slider definition:
sld = javax.swing.JScrollBar(0,1,1,1,size(res,1)+1);		%JScrollBar(int orientation, int value, int extent, int min, int max)
javacomponent(sld, [pos(3)*0.65,5,200,15], gcf);
sld.setUnitIncrement(1); sld.setBlockIncrement(1);
hsld = handle(sld,'CallbackProperties');
set(hsld,'AdjustmentValueChangedCallback',@(source,callbackdata) set(pl,'userdata',res(round(source.Value),:)));
hPl = handle(pl); hProp = findprop(hPl,'UserData');
hPl.addlistener(hProp,'PostSet',@(src,event) refreshdata(gcf));
hPl.addlistener(hProp,'PostSet',@(src,event) title([pl.UserData{2},com,pl.UserData{4}]));
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

%% TDVP time extrapolation
n = 1:para.timeslice; clf;
% n=1:648;
n = n(para.tdvp.calcTime(n,1)>0);
plot(n,para.tdvp.calcTime(n,1)); hold on;
s = csaps(n,para.tdvp.calcTime(n,1));
sn = fnxtr(s);
fnplt(sn,[0,60],0.5,'r')
fnval(sn,60)

%% TDVP compare with Orth2010: 13a -- not really good. Extraction was not the best
% http://arohatgi.info/WebPlotDigitizer/app/?
orthData = cell(5,1);
orthData{1} = load('./../../Orth2010-13a-0.01.csv');
orthData{2} = load('./../../Orth2010-13a-0.05.csv');
orthData{3} = load('./../../Orth2010-13a-0.1.csv');
orthData{4} = load('./../../Orth2010-13a-0.15.csv');
orthData{5} = load('./../../Orth2010-13a-0.2.csv');

myData{1} = res{15,1};
myData{2} = res{16,1};
myData{3} = res{17,1};
myData{4} = res{18,1};
myData{5} = res{19,1};

fignum = 1; clf; hold all;
cellfun(@(x) plot(x(:,1),x(:,2)), orthData, 'UniformOutput', false)
cellfun(@(x) plot(x.para.tdvp.t(:),x.tresults.spin.sz(:)), myData, 'UniformOutput', false)
set(gca,'xlim',[0,325]);
formatPlot(1)

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

%% Format <n> plot for Poster:
ax = gca; f=gcf;
% set(gca, 'LooseInset', [0,0,0,0]); % not with colorbar!!
ax.FontSize = 14;
ax.TickDir = 'out';
% ax.Units = 'centimeters';
% for figure 840x420:
% oldPos = ax.Position;
% ax.Position = [oldPos(1),oldPos(2), 15, 9.8];
f.Units = 'centimeters';
oldPos = f.Position;
f.Position = [oldPos(1),oldPos(2), 19.35, 12.28];	% in cm also produces 15x9.8cm axes
f.Units = 'pixels';
% f.Position = [oldPos(1),oldPos(2), 840, oldPos(4)];
% export_fig('test','-transparent','-png','-m4')
%% inset for poster:
if exist('axIn','var')
	axIn.delete;
end
ax=gca; oldPos = ax.Position;
axIn = axes('Position',[oldPos(1)+0.32, oldPos(2)+0.2, oldPos(3)*0.5, oldPos(4)*0.7]);		% use this if main with colorbar
% axIn = axes('Position',[oldPos(1)+0.36, oldPos(2)+0.2, oldPos(3)*0.5, oldPos(4)*0.7]);	% use this if main without colorbar
axIn.TickDir = 'out';
axIn.FontSize=14; axIn.XColor=[0.7,0.7,0.7];axIn.YColor=axIn.XColor;

hold all;
% mode = 0;		% 0: lin, 1: log
% if mode
% 	surf(tresults.star.omega,tresults.star.t,log10(tresults.star.n));
% 	zlabel('$\log_{10}\left<n_k\right>$');
% else
% 	surf(tresults.star.omega,tresults.star.t,tresults.star.n);
% 	zlabel('$\left<n_k\right>$');
% end
surf(tresults.star.omega,tresults.star.t, ax.Children.ZData);
% xlabel('Mode $\omega_k / \omega_c$');
% ylabel('Time $\omega_c t$');
set(gca,'View',[0 90]);
shading interp
axis tight
axIn.XLim = [0,0.2]; axIn.YLim = [0,420];
[M,I] = max(ax.Children.ZData(:,:),[],2);
plot3(tresults.star.omega(I(2:end)),tresults.star.t(2:end),ones(1,length(tresults.star.t)-1).*10,':r')
plot3([0.1,0.1],[0,350],[10,10],':black');

%% Analyze Chain Parameters SBM_genpara.m

para.s = 1; para.alpha=0.2; para.L=10; para.chain.spectralDensity = 'Leggett_Hard'; para.chain.mapping='OrthogonalPolynomials';para.chain.method='Analytic';
para1 = SBM_genpara(para);
% para.z=1; para.Lambda=1.01; para.chain.mapping = 'OrthogonalPolynomials';para.chain.discrMethod = 'Numerical';para.chain.discretization = 'Linear';para.chain.method = 'Stieltjes';
% para2 = SBM_genpara(para);
%%
figure(1);clf;hold on;
% scatter(para1.epsilon,para1.t); %set(gca,'yscale','log');
plot([para1.epsilon, para1.t]);

%% Map from chain to star:
obs = getObservable({'bath2correlators'},mps,Vmat,para);
fignum = 5;
figure(fignum);clf;
surf(real(obs)); xlabel('Site m'); ylabel('Site n'); set(gca,'view',[0,90]);
title('$Re(\left<a^\dagger_m a_n\right>)$')
fignum = fignum+1; figure(fignum);clf;
surf(imag(obs)); xlabel('Site m'); ylabel('Site n'); set(gca,'view',[0,90]);
title('$Im(\left<a^\dagger_m a_n\right>)$')
fignum = fignum+1; figure(fignum);clf;
surf(abs(obs)); xlabel('Site m'); ylabel('Site n'); set(gca,'view',[0,90]);
title('$\left|\left<a^\dagger_m a_n\right>\right|$')

%%
obs = getObservable({'staroccupation'},mps,Vmat,para);
figure(10);clf;
plot(obs(1,:),obs(2,:));