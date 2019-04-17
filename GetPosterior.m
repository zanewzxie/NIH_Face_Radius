% GetPosterior


% load file
%     load Fitted_Exp5_Pos_Pos_Neu.mat
    load Fitted_Exp4_Neg_Neg_Neu.mat

%%
x = samples.Pmconraw(:,1)-samples.Pmconraw(:,2);
y = samples.SDconraw(:,1)-samples.SDconraw(:,2);
% 
PMpd = fitdist(x,'tLocationScale');
SDpd = fitdist(y,'tLocationScale');
% 
% PMpd = fitdist(x,'Kernel','BandWidth',0.005);
% SDpd = fitdist(y,'Kernel','BandWidth',0.55);

% PMpd = fitdist(x,'Kernel','Kernel','epanechnikov');
% SDpd = fitdist(y,'Kernel','Kernel','epanechnikov');
%%
x_values = -0.2:0.001:0.2;
y_values = -20:0.1:20;

PmPosterior = pdf(PMpd,x_values);figure(3);plot(x_values,PmPosterior.* min(diff(x_values)),'b')
SDPosterior = pdf(SDpd,y_values);figure(4);plot(y_values,SDPosterior.* min(diff(y_values)),'b')

Pm_PosteriorAtNull = PmPosterior(x_values==0).*min(diff(x_values));
SD_PosteriorAtNull = SDPosterior(y_values==0).*min(diff(y_values));
