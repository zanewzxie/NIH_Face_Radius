%clear all; 
close all;
for i =1:10000
    abound = unifrnd(0, 0.35, 1, 1);
    bbound = unifrnd(abound, 0.35, 1, 1);
    sd = gamrnd(0.5,0.5,1,1);
    a= abound; 
    b= bbound;
    
    mu1(i) = normrnd(a, sd, 1, 1);       
    mu2(i) = normrnd(b, sd, 1, 1); 
    mudiff(i)=mu2(i)-mu1(i);
    
    
    % if a < b
    %     mu1(i)=a;
    %     mu2(i)=b;
    % elseif a > b
    %     mu1(i)=b;
    %     mu2(i)=a;
    % else
    %     mu1(i)=NaN;
    %     mu2(i)=NaN;
    % end
    %
    %
    % a = normrnd(0,0.1,1);
    % b = normrnd(0,0.1,1);
    
    priormu1(i)=a;
    priormu2(i)=b;
    ss(i)=sd;
end

diffmu = priormu2-priormu1;
MdiffMu = rad2deg(mean(diffmu))
figure;
subplot(1,3,1); hist(rad2deg(priormu1), 30)
subplot(1,3,2); hist(rad2deg(priormu2), 30)
subplot(1,3,3); hist(rad2deg(diffmu), 30)

figure;
subplot(1,3,1); hist(rad2deg(mu1), 30)
subplot(1,3,2); hist(rad2deg(mu2), 30)
subplot(1,3,3); hist(rad2deg(mudiff), 30)
