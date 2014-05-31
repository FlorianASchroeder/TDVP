function PlotData = ExtractAllData()
% Used to get useful data from results.mat of each subfolder.
% saves extractes data into file:
% Makes Plot to visualize data
%
% Modified:
%      FS 24/05/2014:   Started to include extraction of parameters for MLSBM
%%
% Overwrite extracted.mat ?
overwrite = 0;
global loadDataFromFile
loadDataFromFile = '/results.mat';
global saveDataToFile
saveDataToFile = '/extracted.mat';
global saveDataToFileError
saveDataToFileError = '/extractedError.mat';

regexpLINUX = ['[^:;]*'];
regexpWIN = ['[^;]*'];

doPlot = 0;

if strfind(system_dependent('getos'),'Windows')
    regexpUSE = regexpWIN;
elseif strfind(system_dependent('getos'),'Linux')
    regexpUSE = regexpLINUX;
else
    disp('OS not recognised');
end
dirs = regexp(genpath(pwd),regexpUSE,'match');      % get all subfolders of current directory
dirs = dirs(2:end);                                 % remove current folder

%%
for i = 1:1:length(dirs)
    if (overwrite == 0 && exist([dirs{i},saveDataToFile],'file')) || ~exist([dirs{i},loadDataFromFile],'file')
        continue
    end
    % Extract Data and save to File, per folder
	ExtractOneData(dirs{i});
end

%% Load already existing PlotData??


%% generate Parameter and Essentials list
PlotData = zeros(length(dirs),10);   % (s, alpha, log(shifterror))
PPCWavefunction = zeros(length(dirs),16);   % only applies to Rs.Molischianum
for i = 1:1:length(dirs)
    try
%        load([dirs{i},'/extracted.mat']);
        load([dirs{i},'/results.mat']);         % for MLSBM
        if ~strcmp(para.model,'MLSpinBoson')
            PlotData(i,1) = para.s;
            PlotData(i,2) = para.alpha;
            PlotData(i,3) = -log10(extracted.ShiftMaxError);
            PlotData(i,4) = -log10(extracted.NMaxError);
            PlotData(i,5) = extracted.timeMins;
            PlotData(i,6) = extracted.Sx;       % <sx>
            PlotData(i,7) = para.hx;            % = delta
            PlotData(i,8) = results.spin.sz;    % <sz>
            if max(para.dk) < 1000
                PlotData(i,9) = norm(abs(extracted.aaCreationCorrelator(1:end-1)-extracted.aaCreationCorrelator2));
            else
                PlotData(i,9) = -1;     % dimension too high--> no calculation
            end
            PlotData(i,10) = max(abs(results.nx)); % Max(<n>) to see changes in chain population
        else
            if ~isfield(results,'time');
                continue;
            end
            PlotData(i,1) = para.Lambda;
            PlotData(i,2) = para.z;
            PlotData(i,3) = results.participation;  % extracted.participation; = inverse participation of Cogdell 2006
            PlotData(i,4) = max(abs(results.nx));   % Max(<n>) to see changes in chain population
            PlotData(i,5) = results.E;
            PlotData(i,6) = para.MLSB_p;            % extracted.period;
            if isfield(para,'MLSB_etaFactor')
                PlotData(i,7) = para.MLSB_etaFactor;
            else
                PlotData(i,7) = 1;      % case before MLSB_etaFactor was introduced.
            end
            PlotData(i,8) = para.loop;
            % copy the wavefunction of the ring
            PPCWavefunction(i,:) = diag(calReducedDensity(mps,Vmat,para,1));
        end
    catch

    end

    clear('para','extracted','results');
end
%dlmwrite('20140310-Delta-Alpha-Scan.txt',sortrows(PlotData,7));
save([datestr(now,'yyyymmdd'),'-ExtractedResults.mat'],'PlotData','PPCWavefunction');

if ~doPlot
	return;
end

%% Make Plots if doPlot = 1;

colormap jet
figure(1)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,3),40,PlotData(:,3),'fill')
title('Deviation of Shifts');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$log_{10}(\Delta \delta)$');
figure(2)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,4),40,PlotData(:,4),'fill')
title('Deviation of $<n>$');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$log_{10}(\Delta <n>)$');
figure(3)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,5),40,PlotData(:,5),'fill')
title('Time needed');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('t in min');
figure(4)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,6),30,PlotData(:,7),'fill')
title('$<\sigma_x>$ vs $\alpha$');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$<\sigma_x>$');
%% Plot <sz> vs delta, alpha
figure(5)
subplot(1,2,1)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,8),10,PlotData(:,8),'fill')
title('$<\sigma_z>$ vs $\alpha$');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$<\sigma_z>$');
subplot(1,2,2)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,8),10,PlotData(:,8),'fill')
set(gca,'View',[0 90]);
title('$<\sigma_z>$ vs $\alpha$');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$<\sigma_z>$');
%% Plot error in <a_i a_i+1> vs delta, alpha
figure(6)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,9),10,PlotData(:,9),'fill')
title('$<a_i^\dagger a_{i+1}^\dagger> - (f_i f_{i+1}-<a_i^\dagger><a_{i+1}^\dagger>)$ vs $\alpha$');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$<\sigma_z>$');
%% Plot Max(<n>) vs delta, alpha
figure(7)
scatter3(PlotData(:,7),PlotData(:,2),PlotData(:,10),10,PlotData(:,10),'fill')
title('$Max(<n>)$ vs $\alpha$');
xlabel('$\Delta$');
ylabel('$\alpha$');
zlabel('$Max(<n>)$');
%% MLSBM from here:
%  PlotData: 1 = Lambda; 2 = z; 3 = participation; 4 = Max(<n>); 5 = E0; 6 = period;
%            7 = etaFactor;

figure(1)
scatter3(PlotData(:,1),PlotData(:,2),PlotData(:,3),10,PlotData(:,3),'fill')
title('$Max(<n>)$ vs $\Lambda$ and z');
xlabel('$\Lambda$');
ylabel('$z$');
zlabel('$Max(<n>)$');
%% Plot Participation for entire Dataset.
figure(2)
plot(PlotData(:,3));
%% Plot Max Chain Occupation vs. Period and Strength
figure(2)
scatter3(PlotData(:,6),PlotData(:,7),PlotData(:,4),50,PlotData(:,4),'fill')
xlabel('$p$');
ylabel('$\eta$')
zlabel('$<n_1>$');
%% Plot Loop vs. Period and Strength
figure(2)
scatter3(PlotData(:,6),PlotData(:,7),PlotData(:,8),50,PlotData(:,8),'fill')
xlabel('$p$');
ylabel('$\eta$')
zlabel('Loop');
formatPlot(2)
%% Plot Participation vs. Period and Strength
figure(2)
scatter3(PlotData(:,6),PlotData(:,7),PlotData(:,3),50,PlotData(:,3),'fill')
xlabel('$p$');
ylabel('$\eta$')
zlabel('$P$');
formatPlot(2)
%% Plot Participation vs. Period
figure(2)
%select = PlotData(:,6)==4;          % period == 4
select = PlotData(:,7)==1;          % etaFactor == 4
plot(PlotData(select,6),PlotData(select,3))
xlabel('$p$');
ylabel('$P$');
%% Plot Participation vs. Strength
defPeriods = [2,3,4,5,7,8,9,15,16];
colors=jet(length(defPeriods));
figure(2)
hold all
i=1;
for periods = defPeriods
    select = PlotData(:,6)==periods & PlotData(:,7)>=0;          % period == 4 or 16
    Indexed = [PlotData(select,7) PlotData(select,3)];
    Indexed = sortrows(Indexed,1);
    plot(Indexed(:,1),Indexed(:,2),'Marker','*','LineStyle','none','Color',colors(i,:));
    i=i+1;
end
xlabel('$\eta$');
ylabel('$P$');
legend({'2','3','4','5','7','8','9','15','16'})
formatPlot(2);
line([0 5],[16 16],'LineWidth',1,'Color','black');
%% Plot GS Energy vs. Period and Strength
figure(3)
scatter3(PlotData(:,6),PlotData(:,7),PlotData(:,5),50,PlotData(:,5),'fill')
xlabel('$p$');
ylabel('$\eta$')
zlabel('$E_0$');

%% Plot the wavefunction for different strengths
% modify the required period
select =  PlotData(:,6)==16;
IndexedMatrix = [PlotData(select,7) PPCWavefunction(select,:)];
IndexedMatrix = sortrows(IndexedMatrix,1);                      % sort with respect to y axis
figure(4);
surf(1:16,IndexedMatrix(:,1),IndexedMatrix(:,2:end));
xlabel('Site $k$');
ylabel('$\eta$')
zlabel('$|\Psi(k)|^2$')
formatPlot(4)

%% Plot the wavefunction for different periods
select = PlotData(:,7)==0; %& PlotData(:,7) >=0;
IndexedMatrix = [PlotData(select,6) PPCWavefunction(select,:)];
IndexedMatrix = sortrows(IndexedMatrix,1);
figure(4);
surf(1:16,IndexedMatrix(:,1),IndexedMatrix(:,2:end));
xlabel('Site $k$');
ylabel('$p$')
zlabel('$|\Psi(k)|^2$')

%% Scatter interpolation
%% using griddata
figure
[xi, yi] = meshgrid(1:0.5:29, -1:0.1:5);
x = PlotData(:,6);
y = PlotData(:,7);
z = PlotData(:,3);
zi = griddata(x,y,z, xi,yi);
surf(xi,yi,zi);
xlabel('$p$'), ylabel('$\eta$'), zlabel('$P$')
%% using Interpolant class
F = scatteredInterpolant([x y],z)
[Xq,Yq] = meshgrid(1:0.5:40, -1:0.1:5);
F.Method = 'natural';
Vq = F(Xq,Yq);
surf(Xq,Yq,Vq);
xlabel('X','fontweight','b'), ylabel('Y','fontweight','b');
zlabel('Value - V','fontweight','b');
title('Linear Interpolation Method','fontweight','b');
end

function ExtractOneData(path)
% Extracts Data and saves it to file
% Use as:
% ExtractOneData(dirs{i});
global loadDataFromFile
global saveDataToFile
global saveDataToFileError
try
	load([path,loadDataFromFile]);
catch
	fprintf('%s : file does not exist\n',path);
	return;
end
E = struct();
if isfield(results,'time');
	E.timeMins = results.time / 60;
else
	E.aborted = 1;
end
E.loops = para.loop;
%% Calculate Shift analytically for iSBM:
    t = para.t; e = para.epsilon;
    A = gallery('tridiag',t(2:end),e(1:end),t(2:end));      % creates tridiag for system A.(shiftVec) = (-t1*sigmaZ, 0,0,0,0,...)
    B = zeros(para.L-1,1);
    if isfield(results,'spin')                              % if spin was calculated in the end
        B(1) = -t(1)*(results.spin.sz);
    else
        B(1) = -t(1)*1;
    end
    % E.ShiftCalc = results.shift{1};                        % this is analytical solution for spin.
    E.ShiftCalc = [0; A\B.*sqrt(2)]';
%% Compare analytical and simulated iSBM shift
    E.ShiftSim = para.shift;
    E.ShiftDeviation = E.ShiftCalc./E.ShiftSim-1;
    E.ShiftMaxError = max(abs(E.ShiftDeviation));           % Max absolute Error
    E.ShiftMeanError = mean(abs(E.ShiftDeviation(2:end)));  % Mean absolute Error
    E.ShiftMax = max(abs(E.ShiftCalc));                     % Max Shift
%% Compare Occupation numbers
    E.NCalc = (E.ShiftCalc.*E.ShiftCalc./2);
    try
        E.NSim = results.nx;
        E.converged = 1;
    catch
        E.converged = 0;
        fprintf('%s is not converged\n',para.filename(1:58));
        return;
        if max(para.dk) < 1000
            E.NSim = calbosonocc_SBM1(mps,Vmat,para,results);   % still calculate from approximate state, takes quite long in bad cases!
        else
            fprintf('%s has dk = %u\n',para.filename(1:58),max(para.dk));
            return;
        end
    end
    E.NDeviation = E.NCalc./E.NSim-1;
    E.NMaxError = max(abs(E.NDeviation));
    E.NMeanError = mean(abs(E.NDeviation(2:end)));
    E.NMaxCalc = max(abs(E.NCalc));
    E.NMaxSim = max(abs(E.NSim));
%% Save spin <sx> for Plot against alpha
% try because not applicable to MLSBM
    if E.converged
        try
            E.Sx = results.spin.sx;
        catch
        end
    else
        try
            E.spin = calspin(mps,Vmat,para,results);
            E.Sx = E.spin.sx;
        catch
        end
    end
%% Calculate <a_i a_i+1>, call aaCreationCorrelator
    if max(para.dk) < 1000
        E.aaCreationCorrelator = calboson2siteCreationCorrelator_SBM1(mps,Vmat,para);   % calculated directly
        E.aCreationCorrelator = calboson1siteCreationCorrelator_SBM1(mps,Vmat,para);    % <a^+_i>
        E.aaCreationCorrelator2 = para.shift(1:end-1).*para.shift(2:end) - E.aCreationCorrelator(1:end-1).*E.aCreationCorrelator(2:end);            % calculated as f_i f_i+1 - <a^+_i><a^+_i+1>
    else
        fprintf('%s has dk = %u\n',para.filename(1:58),max(para.dk));
        return;
    end

%% Results specifically for MLSBM
    if strcmp(para.model,'MLSpinBoson')
        try
            E.participation = results.participation;
        catch
            E.participation = calParticipation(calReducedDensity(mps,Vmat,para,1));
        end
		E.EvaluesLog = cell2mat(results.EvaluesLog);
		E.EnergyConvergence = E.EvaluesLog-min(E.EvaluesLog);					% Plot against log-axis therefore -min
        E.E0 = results.E;
        E.period = para.MLSB_p;
    end

%% Save to file
extracted = E;
if E.converged
	save([path,saveDataToFile],'extracted','para','results');
else
	save([path,saveDataToFileError],'extracted','para','results');
end
end
