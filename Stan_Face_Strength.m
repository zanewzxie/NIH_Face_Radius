
% add stan toolbox. 
addpath(genpath('/Users/macpro/Documents/MATLAB/MatlabStan/'), '-end')
addpath(genpath('/Users/macpro/Documents/Data/Hyungbum/BHM_Park_Play/MemToolbox'), '-end')


%% Fit mulitple subject's data using Stan for various conditions All at once!
clear all; close all

if exist('temp.data.R', 'file') ~=0
    delete('temp.data.R')
end

if exist('output-1.csv', 'file') ~=0
    delete('output-1.csv')
end
if exist('output-2.csv', 'file') ~=0
    delete('output-2.csv')
end
if exist('output-3.csv', 'file') ~=0
    delete('output-3.csv')
end
if exist('output-4.csv', 'file') ~=0
    delete('output-4.csv')
end

% import data
load('/Users/Zane/Documents/Bethesda/Research/NIH_Face_Radius/Exp4Face.mat')
%load('TwoVM_Exp3.mat')
%% Simulation

% clear all
% model.name = '2VMu';
% model.paramNames = {'mu1','mu2','sd1', 'g'};
% model.lowerbound = [-180 -180 0 0];
% model.upperbound = [180 180 360 1];
% model.movestd = [0.05, 0.05, 0.05, 0.01];
% model.pdf = @(data, mu1, mu2, sd1, g) ((1-g).*(0.5.*vonmisespdf(data.errors(:),mu1,deg2k(sd1)) + 0.5.*vonmisespdf(data.errors(:),mu2,deg2k(sd1))) + ...
%     g.*1/360);
% model.start = [0 10 20 0.2; 10 20 25 0.15; -10 5 20 0.1; -5 5 20 0.3];
% 
% paramset = [10 15 25 0.2];
% NumSimSub = 2; NumSimTrial = 500;
% for i = 1:NumSimSub
%     SimData((i-1)*NumSimTrial + 1: i*NumSimTrial, 1) = i; 
%     SimData((i-1)*NumSimTrial + 1: i*NumSimTrial, 2) = SampleFromModel(model, paramset, [NumSimTrial, 1]);
% end
% TwoVM = SimData; %% temp

%% Make data structure as Data.errors(ncond, nsub, ntrial)
nsub =  10; % length(unique(TwoVM(:,1)));
% ntrials = length(TwoVM)/nsub;
Data.nsub  = nsub;
Data.ncond  = 4; 
Data.RadiusLevel = [2 3 4 5];

Data.allerrors=[];
Data.allsubj=[];
Data.allradius=[];

for isub=1:nsub
    mask=DataAll{isub}.Task==1 & DataAll{isub}.WhichWheel>=2;   
    Data.allerrors=[ Data.allerrors; DataAll{isub}.errors(mask)];
    Data.allradius=[ Data.allradius; DataAll{isub}.WhichWheel(mask)];
    Data.allsubj=[ Data.allsubj; ones(length(DataAll{isub}.errors(mask)),1)*isub];
end
Data.ntrials=length(Data.allerrors);
a = Data.allerrors;
a = deg2rad(a);
a(a > 3.14159) = 3.14159;
a(a < -3.14159) = -3.14159;
Data.allerrors=a;

%% Stan

% model = 'RAC_Beta_Gamma_OneVM.stan';
% model = 'RAC_Beta_Gamma_TwoVM.stan';
% model = 'RAC_TwoVM_NarrowingTwoMu.stan';
model = 'RAC_TwoVM_OrderedTwoMuUnif.stan';
% model = 'RAC_Beta_Gamma_TwoVM_Full.stan';
niteration = 2000; % for publication purpose;

fit = stan('file',model,'data',Data,'file_overwrite',true,'verbose',true,'iter',niteration,'warmup', 2000, 'refresh', round(max(niteration/20,1)));

%%

%%% 4.5 hrs for 2000-3000 for 50 Subj / 500t
% fit.stop();

print(fit)
samples = fit.extract('permuted',true);
% figure(100);
%%
% close all
% PlotVariable = -1 * rad2deg(samples.GrPdiff_mu);
PlotVariable = [rad2deg(samples.Posterior_GrpMuDiff), rad2deg(samples.Prior_GrpMuDiff)];
% PlotVariable = [-1*rad2deg(samples.GrPdiff_mu), -2*rad2deg(samples.Prior_GrPdiff_mu)];
FigureLabel = {'Posterior', 'Prior'};

h=[1;1]; FigHandle = figure; set(FigHandle, 'Position', [50 600 700 500])
ncond = size(PlotVariable,2);
color ={'r','b','g'};

%Overlapping distributions
subplot(2,1,1);
hist([rad2deg(samples.GrpMu1), rad2deg(samples.GrpMu2)], 30); 
title('Posterior Distribution of GrpMu1 & GrpMu2','fontsize',16);
legend({' GrpMu1', ' GrpMu2'}, 'box','off','fontsize',12, 'Location','northeast'); xlim([5, 20]);

subplot(2,1,2)
MuX = [-5:0.01:30];%rad2deg(0.35)]; 
PosDist = fitdist(PlotVariable(:,1), 'Kernel','BandWidth', 0.5); %'Kernel','BandWidth', 0.7); %'tlocationscale');	%'tlocationscale');
PriorDist = fitdist(PlotVariable(:,2), 'half normal'); %'Kernel','BandWidth',0.9); 

    %PriorDist = unifpdf(MuX, 0, 10);
    %PriorDist = fitdist(zzz, 'Kernel','BandWidth', 0.5); 

% MuDiffPos = smooth(pdf(PosDist, MuX),1000); 
MuDiffPos = pdf(PosDist, MuX);
% plot(MuX, MuDiffPos.* min(diff(MuX)),'b', 'LineWidth',2); hold
plot(MuX, MuDiffPos, 'Color', [0.75 0 0], 'LineWidth',1.5); hold on;

MuDiffPrior = pdf(PriorDist,MuX); 
% plot(MuX,MuDiffPrior.* min(diff(MuX)),'k--', 'LineWidth',2)
plot(MuX,MuDiffPrior,'-', 'Color', [0 0 0], 'LineWidth',1.5)

% Median Posterior
HDI = hdiFromSamples(PlotVariable(:,1), 0.95)
H_HDI = HDI(1);  L_HDI=HDI(2);

X = MuX(MuX>=H_HDI & MuX<= L_HDI);
Y = MuDiffPos(MuX>=H_HDI & MuX<= L_HDI);

h(2,:) = area(X,Y,'FaceColor', [240 225 225]/255 ,'LineStyle',':'); 
% plot(median(PlotVariable(:,1)), 0,'ro','markersize',8,'MarkerFaceColor',[0.8 0 0]);
plot(mean(PlotVariable(:,1)), 0,'ro','markersize',8,'MarkerFaceColor',[0.8 0 0]);

plot(MuX, MuDiffPos, 'Color', [0.75 0 0], 'LineWidth',1.5);
plot(MuX,MuDiffPrior,'-', 'Color', [0 0 0], 'LineWidth',1.5)
%[f1,x1] = ksdensity(PlotVariable(:,2)-PlotVariable(:,1),'NumPoints' ,3000,'bandwidth',0.5,'kernel','epanechnikov');

MuDiffPosNull = MuDiffPos(MuX==0).*min(diff(MuX));
MuDiffPriorNull = MuDiffPrior(MuX==0).*min(diff(MuX));

title('Posterior Distribution of GrpMu1 & GrpMu2','fontsize',16);
legend({'cPosterior GrpDiffMu', ' Prior GrpDiffMu', ' 95% HDI', ' Median_G_r_p_D_i_f_f_?'}, 'box','off','fontsize',12, 'Location','northeast'); xlim([-5, 30]);

% figure;
% subplot(2,3,1); hist(rad2deg(samples.Prior_GrpMu1),15)
% subplot(2,3,2); hist(rad2deg(samples.Prior_GrpMu2),15)
% subplot(2,3,3); hist(rad2deg(samples.Prior_GrpMuDiff),15)
% 
% subplot(2,3,4); hist(rad2deg(samples.GrpMu1),15)
% subplot(2,3,5); hist(rad2deg(samples.GrpMu2),15)
% subplot(2,3,6); hist(rad2deg(samples.Posterior_GrpMuDiff),15)

BF10 = MuDiffPriorNull/MuDiffPosNull
% save Results_2vmFull1SD.mat fit samples Data model






set(gca,'fontsize',12);
ylabel('Density','FontSize',13); xlabel('Shift (Mu)','FontSize',13);
title('Posterior Distribution of Mu Difference','fontsize',16);
legend(FigureLabel, 'box','off','fontsize',14, 'Location','northeast')
%axis([min(min(x)) 0 0 ceil(max(max(f)))])
%axis AUTO%([min(MuX) ceil(max(MuX)) 0 0 ])

mxlp=find(samples.lp__==max(samples.lp__));
MAP.PmSub = samples.pm(mxlp,:); 
MAP.Mu1Sub = rad2deg(samples.mu1(mxlp,:)); MAP.Mu2Sub = rad2deg(samples.mu2(mxlp,:)); MAP.SDSub = rad2deg(k2sd(samples.kappa1(mxlp,:)));
MAP.Mu1 = rad2deg(samples.GrpMu1(mxlp,:)); MAP.Mu2 = rad2deg(samples.GrpMu2(mxlp,:)); 
MAP.MuDiff = rad2deg(samples.Posterior_GrpMuDiff(mxlp,:));
MAP.SD = mean(rad2deg(k2sd(samples.kappa1(mxlp,:))));
MAP.Pm = mean(samples.pm(mxlp,:));

Mparam.Mu1 = mean(rad2deg(samples.GrpMu1)); Mparam.Mu2 = mean(rad2deg(samples.GrpMu2)); Mparam.MuDiff = mean(rad2deg(samples.Posterior_GrpMuDiff));
Mparam.Pm = mean(mean(samples.pm)); Mparam.SD = mean(mean(rad2deg(k2sd(samples.kappa1))))


% subplot(1,3,3); hist(SimData(:,2),30);
% save Results_temp.mat fit samples Data model

%% More analysis and ploting

% 1. Plot the posterior distribution with HDI and mode mean being marked.

%clear all; close all; load Results_Temp
addpath(genpath('/Applications/ToolBox/BHM_Park_Play/pmtk3'));
% Calculate 95% 
FigureLabel = {'GrpMu1', 'GrpMu2'};


figure(1); 
subplot(1,2,1); hist([rad2deg(samples.GrpMu1), rad2deg(samples.GrpMu2)], 30); xlim([0, 20]);
subplot(1,2,2); histfit(rad2deg(samples.Posterior_GrpMuDiff),20,'kernel'); xlim([0, 10])
set(figure(1), 'Position', [50 550 1000 250]); 

% PlotVariable = -1 * rad2deg(samples.GrPdiff_mu);
PlotVariable = [rad2deg(samples.GrpMu1), rad2deg(samples.GrpMu2)];

h=[1;1]; FigHandle = figure; set(FigHandle, 'Position', [50 200 1000 250])
ncond = size(PlotVariable,2);
color ={'r','b','g'};

%Overlapping distributions
subplot(1,2,1)
for i = 1:ncond
    [f(:,i),x(:,i)] = ksdensity(PlotVariable(:,i), 'NumPoints',3000,'bandwidth',0.5);
    h(:,i)= plot(x(:,i),f(:,i),[color{i} '-'],'LineWidth',2); hold on
    %axis([min(x(:,i))-0.15 max(x(:,i))+0.1 0 20]);set(gca, 'YTick', []);
end
set(gca,'fontsize',12); xlim([0, 20]); ylim([0, 0.6]);
ylabel('Density','FontSize',13); xlabel('Shift (Mu)','FontSize',13);
title('Posterior Distribution of Mu','fontsize',13);
legend(h(1,:), FigureLabel,'box','off','fontsize',14, 'Location','northwest')


subplot(1,2,2)
%HDI = hdiFromSamples(PlotVariable(:,2)-PlotVariable(:,1),0.95);
HDI = rad2deg(hdiFromSamples((samples.Posterior_GrpMuDiff), 0.95))
H_HDI = HDI(1);  L_HDI=HDI(2) ;
[f1,x1] = ksdensity(PlotVariable(:,2)-PlotVariable(:,1),'NumPoints' ,3000,'bandwidth',0.5,'kernel','epanechnikov');
    %[f2,x2] = ksdensity(c,'NumPoints' ,3000,'bandwidth',5,'kernel','epanechnikov');
X = x1(x1>=H_HDI & x1<= L_HDI);
Y = f1(x1>=H_HDI & x1<= L_HDI);
h(2,:) = area(X,Y,'FaceColor', [0.8 0.95 0.95],'LineStyle',':')  ;hold on;
h(2,:)= plot(x1,f1,'k-','LineWidth',2); axis([min(x1)-0.1 max(x1)+0.1 0 max(f1)+1]); hold on;
    %plot(x2,f2,'k-','LineWidth',2); axis([min(x2)-0.1 max(x2)+0.1 0 max(f2)+1]); hold on;
plot(mean(PlotVariable(:,2)-PlotVariable(:,1)),0,'ro','markersize',8,'MarkerFaceColor',[0.75 0 0])

set(gca,'fontsize',12); xlim([-5, 10]); ylim([0, 0.6]); % set(gca, 'YTick', []); 
ylabel('Density','FontSize',13); xlabel('Shift (Mu)','FontSize',13);
title('Posterior Distribution of MuDiff','fontsize',13);
legend(h(1,:), FigureLabel,'box','on','fontsize',16, 'Location','northwest')


%'FontWeight','bold'
% export_fig Affect_LTM_Exp3_Pm -pdf

