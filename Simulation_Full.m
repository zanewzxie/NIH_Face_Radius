% clear all;
% 
% load Results_2vmFull2SD.mat
cd('/Users/macpro/Documents/Data/Hyungbum/BHM_Park_Play/Consolidation/Results/NewResults')
niteration = size(samples.pm, 1);
nsub = size(samples.pm, 2);
ntrial = size(Data.errors,1);

%% 2VM 2VM 2VM
Break = 2;
for i = 1:niteration/Break
    Remaining = 100-((niteration-i)/niteration *100)
    for isub = 1:nsub
        pm(i,isub)  = samples.pm(i,isub);
        kappa1(i,isub) = samples.kappa1(i,isub);
%         kappa2(i,isub) = samples.kappa2(i,isub);
        mu1(i,isub) = samples.mu1(i,isub);
        mu2(i,isub) = samples.mu2(i,isub);
        SimData = Data.errors(:,isub);
        
        Simmodel.pdf = @(data, mu1, mu2, kappa1, kappa2, pm) (pm.*(0.5.*vonmisespdf(data(:),mu1,kappa1) + 0.5.*vonmisespdf(data(:),mu2,kappa2)) + (1-pm).*(1/360));
        
        for itrial = 1:ntrial
            Likelihood(i,isub,itrial) = Simmodel.pdf(SimData(itrial), mu1(i,isub), mu2(i,isub), kappa1(i,isub), kappa1(i,isub), pm(i,isub));
            loglikelihood(i,isub,itrial) = log(Likelihood(i,isub,itrial));
        end
        
        
    end
    log_lik(i,:) = reshape(loglikelihood(i,:,:),nsub*ntrial,1);
    clc
end

for i = niteration/Break+1:niteration
    Remaining = 100-((niteration-i)/niteration *100)
    for isub = 1:nsub
        pm(i,isub)  = samples.pm(i,isub);
        kappa1(i,isub) = samples.kappa1(i,isub);
%         kappa2(i,isub) = samples.kappa2(i,isub);
        mu1(i,isub) = samples.mu1(i,isub);
        mu2(i,isub) = samples.mu2(i,isub);
        SimData = Data.errors(:,isub);
        
        Simmodel.pdf = @(data, mu1, mu2, kappa1, kappa2, pm) (pm.*(0.5.*vonmisespdf(data(:),mu1,kappa1) + 0.5.*vonmisespdf(data(:),mu2,kappa2)) + (1-pm).*(1/360));
        
        for itrial = 1:ntrial
            Likelihood(i,isub,itrial) = Simmodel.pdf(SimData(itrial), mu1(i,isub), mu2(i,isub), kappa1(i,isub), kappa1(i,isub), pm(i,isub));
            loglikelihood(i,isub,itrial) = log(Likelihood(i,isub,itrial));
        end
        
        
    end
    log_lik(i,:) = reshape(loglikelihood(i,:,:),nsub*ntrial,1);
    clc
end

% log_lik = [log_lik1; log_lik2];




%% 1VM 1VM 1VM
for i = 1:niteration
    
    Remaining = 100-((niteration-i)/niteration *100)
    for isub = 1:nsub
        pm(i,isub)  = samples.pm(i,isub);
        kappa1(i,isub) = samples.kappa(i,isub);
        mu(i,isub) = samples.mu(i,isub);        
        SimData = Data.errors(:,isub);
        
        Simmodel.pdf = @(data, mu, kappa, pm) (pm.*vonmisespdf(data(:),mu,kappa) + (1-pm).*(1/360));
        
        for itrial = 1:ntrial
            Likelihood(i,isub,itrial) = Simmodel.pdf(SimData(itrial), mu(i,isub), kappa1(i,isub), pm(i,isub));
            loglikelihood(i,isub,itrial) = log(Likelihood(i,isub,itrial));
        end
        
        
    end
    log_lik(i,:) = reshape(loglikelihood(i,:,:),nsub*ntrial,1);
    clc
end
%%

%figure; hist(rad2deg(MU_mu), 30)

% WAIC_Exp5_SS1_2VMnew = mstan.waic(log_lik)
WAIC_Exp1_2VM_Mu1Diff = mstan.waic(log_lik)
% SimMuGap4_2VM_M7 = mstan.waic(log_lik)

% save Results_Exp1_SS2_2VMorderNew2.mat TwoVM samples paramset SimMuGap4_2VM_M7
save Results_Exp1_2VM_Mu1Diffwide.mat code samples Data WAIC_Exp1_2VM_Mu1Diff %model fit 
% save Test_2VMorderFull.mat samples Data WAIC_Exp2_2VMorderNew %model fit 
