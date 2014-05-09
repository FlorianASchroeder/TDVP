%% Path where to save figures
path = pwd;
saveto = '..\Presentations\20140206 - Summary\';
wantSave = 0;
%% Plot Vmat contributions for different Sites i (normalized)
i=2;
plot(Vmat{1,i}(:,:))
title(['k = ',num2str(i),', max SV = ',num2str(results.Vmat_sv{1,i}(1,1))])
ylabel('Contribution to OBB')
xlabel('$d_k$')
if wantSave
    export_fig(sprintf('%sVmatNormalized%s-%u',saveto,para.filename(1:13),i),'-transparent','-png','-eps')
end
%% Plot Vmat contributions for different Sites i (SV-weighted)
i = 2;
plot(Vmat{1,i}*diag(results.Vmat_sv{1,i}))
title(['k = ',num2str(i),', max SV = ',num2str(results.Vmat_sv{1,i}(1,1))])
ylabel('Contribution to OBB')
xlabel('$d_k$')
%print(gcf, [saveto,'VmatScaled',num2str(i),'.eps'],'-deps')
if wantSave
    export_fig(sprintf('%sVmatScaled%s-%u',saveto,para.filename(1:13),i),'-transparent','-png','-eps')
end

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
    plot(results.nx);
    title('$$<n_x(k)>$$');
subplot(2,2,2);
    plot(para.trustsite);
    title('Trustsite')
% 3D-version elucidating change:
% subplot(2,2,3);
%     surf(cell2mat(results.shift'))
%     set(gca,'View',[-25 10]);
%     shading interp
%     title('Bosonic shift');
% subplot(2,2,4);
%     surf(cell2mat(results.d_opt'));
%     shading interp
%     set(gca,'View',[0 90]);
%     title('OBB dim')
% 2D-version showing final
subplot(2,2,3);
    plot(results.shift{end});
    title('Bosonic shift');
subplot(2,2,4);
    plot(results.d_opt{end});
    title('OBB dim')
    text(-80,-30,sprintf(para.filename(1:38)))
if wantSave
    export_fig(sprintf('%sResultsSummary%s',saveto,para.filename(1:13)),'-transparent','-png')
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
set(pl(1), 'Marker', 'none','Color',colors{1}); % blue
pl(2) = plot(para.shift);
set(pl(2), 'Marker', 'none','Color',colors{2}); % green
if wantSave
    export_fig(sprintf('%sAnalyticShift%s',saveto,para.filename(1:13)),'-transparent','-png','-eps')
end
%% Plot deviation from calculated occupation
nx = [0 (shift.*shift./2)'];
%plot(nx)
%plot(results.nx)
plot((nx-results.nx));
title('Absolute Deviation of $<n_k>$');
xlabel('site k');
ylabel('$\Delta <n_k>$')
if wantSave
    export_fig(sprintf('%sAbsoluteDeviationN%s',saveto,para.filename(1:13)),'-transparent','-png','-eps')
end
%% Plot relative deviation of shift
relShift = ((shift-para.shift(2:end)')./shift);
pl(1) = plot([0; relShift])
title('Relative deviation of shift from calculation')
xlabel('site k')
ylabel('$\frac{\delta_k-\delta_{k,VMPS}}{\delta_k}$')
if wantSave
    export_fig(sprintf('%sRelativeDeviationShift%s',saveto,para.filename(1:13)),'-transparent','-png','-eps')
end

%% Plot s-alpha relation
figure(5)
x = 0:0.01:1;
y = 0.005;
z = y.^(x./(1-x));
plot(x,z);

PlotData(:,PlotData(:,7)==0.001)