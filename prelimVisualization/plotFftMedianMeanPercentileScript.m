% Once the spectra has been generated using
% wiimoteRecodingsAnalysisScript.m

% WDTid 32
% figure(10);plot(median(spectra(:,1:240),2),'r');hold on;
% figure(10);plot(median(spectra(:,241:end),2),'b');hold off;

% figure(11);plot(mean(spectra(:,1:240),2),'r');hold on;
% figure(11);plot(mean(spectra(:,241:end),2),'b');hold off;

% figure(12);plot(prctile(spectra(:,1:240),95,2),'r');hold on;
% figure(12);plot(prctile(spectra(:,241:end),95,2),'b');hold off;

figure(10);
plot(prctile(spectra(:,1:240),95,2),'r');hold on;
plot(prctile(spectra(:,241:end),95,2),'b');
plot(median(spectra(:,1:240),2),'r','LineWidth',2);
plot(median(spectra(:,241:end),2),'b','LineWidth',2);hold off;

figure(11);
nPlot = 1;
for freqComp = [23,43,127];
    [ftshist, binpos] = hist(spectra(freqComp,1:240));
    subplot(1,3,nPlot); 
    bar(binpos,ftshist,'FaceColor',[1 0 0],'FaceAlpha',.5);hold on;
    [ftshist, binpos] = hist(spectra(freqComp,241:end));
    bar(binpos,ftshist,'FaceColor',[0 0 1],'FaceAlpha',.5);hold off;
    nPlot = nPlot+1;
end