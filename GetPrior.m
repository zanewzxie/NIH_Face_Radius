% Get the prior without data
%

nconds = 2;
nsub = 20;
bound = 2;
niter = 100000;

for iter = 1:niter
    
    SdPmoffset = unifrnd(0.001,2);
    SdSDoffset = unifrnd(0.001,2);
    SdPmsub = unifrnd(0.001,2);
    SdSDsub = unifrnd(0.001,2);
    SdPmbase = unifrnd(0.001,2);
    SdSDbase = unifrnd(0.001,2);
    MeanPmbase = unifrnd(-2,2);
    MeanSDbase = unifrnd(-2,2);
    
    iter;
    for isub = 1:nsub
        
        Pmsub(isub,iter) = normrnd(0,SdPmsub);
        SDsub(isub,iter) = normrnd(0,SdSDsub);
        
        for icon = 1:nconds
            
            
            Pmcon(icon,iter) = normrnd(MeanPmbase,SdPmbase);
            SDcon(icon,iter) = normrnd(MeanSDbase,SdSDbase);
            
            Pmoffset(icon,isub,iter) = normrnd(0,SdPmoffset);
            SDoffset(icon,isub,iter) = normrnd(0,SdSDoffset);
            
            pm(icon,isub,iter) =  inv_logit(Pmcon(icon,iter) + Pmsub(isub,iter) + Pmoffset(icon,isub,iter));
            SD(icon,isub,iter) = 5+(100-5)*inv_logit(SDcon(icon,iter) + SDsub(isub,iter) + SDoffset(icon,isub,iter));
            
        end
    end
    
    for icon = 1:nconds
        Pmconraw(icon,iter) = mean(pm(icon,:,iter));
        SDconrawSD(icon,iter) = mean(SD(icon,:,iter));
    end
    
end

save GetPrior.mat Pmconraw SDconrawSD pm SD
% get pm and sd percondition
%%
load GetPrior.mat

%%
% save GetPrior.mat
x = Pmconraw(1,:)'-Pmconraw(2,:)';
y = SDconrawSD(1,:)'-SDconrawSD(2,:)';

PMpd = fitdist(x,'Kernel','Kernel','epanechnikov');
SDpd = fitdist(y,'Kernel','Kernel','epanechnikov');

% PMpd = fitdist(x,'Kernel','BandWidth',0.01);
% SDpd = fitdist(y,'Kernel','BandWidth',1);
%%
x_values = -0.5:0.001:0.5;
y_values = -40:0.1:40;

PmPrior = pdf(PMpd,x_values);hold on;figure(3);plot(x_values,PmPrior.* min(diff(x_values)),'r')
SDPrior = pdf(SDpd,y_values);hold on;figure(4);plot(y_values,SDPrior.* min(diff(y_values)),'r')

Pm_PriorAtNull = PmPrior(x_values==0).*min(diff(x_values));
SD_PriorAtNull = SDPrior(y_values==0).*min(diff(y_values));


Pm_BF = Pm_PriorAtNull/Pm_PosteriorAtNull;
SD_BF = SD_PriorAtNull/SD_PosteriorAtNull;



