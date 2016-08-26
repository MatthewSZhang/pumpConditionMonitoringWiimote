% Set save figure options (Will probably need to modify path below)
saveFigureOption = false;
if (saveFigureOption)
    % Save plot
    folderName = sprintf('caseStudyPumpWDTid%dAcc%s',WDTid,accelerometerSignal);
    pathName = fullfile('C:\Users\engs1602\research\meetings\smallGroup\20160825ManandharAnalysisCodeUpdate\plots',folderName);
end

% Once the spectra has been generated using
% wiimoteRecodingsAnalysisScript.m
prcVal = 95;

figure(10);
medAb = median(spectra(:,1:240),2);
medNo = median(spectra(:,241:end),2);
prctileAb = prctile(spectra(:,1:240),prcVal,2);
prctileNo = prctile(spectra(:,241:end),prcVal,2);
plot(medAb,'r','LineWidth',2);hold on;
plot(medNo,'b','LineWidth',2);
plot(prctileAb,'r');hold on;
plot(prctileNo,'b');hold off;
grid on;
legend('Median Ab','Median No','95th per Ab','95th per No');
xlabel('Hz');
ylabel('|FFT|');
yMax = max(cat(1,medAb,medNo,prctileAb,prctileNo));
axis([0 size(medAb,1)+5 0 yMax+1]);

if (saveFigureOption)
    % Save plot    
    plotName = fullfile(pathName,sprintf('medianVs95PercentileNoVsAbWdt%dAcc%s',WDTid,accelerometerSignal));
    fig = gcf;
    fig.PaperPositionMode = 'auto';    
    print(plotName,'-dpng','-r0');
    plotName = fullfile(pathName,sprintf('medianVs95PercentileNoVsAbWdt%dAcc%s.fig',WDTid,accelerometerSignal));
    savefig(plotName);
end

% Using hist to approx pdf
% figure(11);
% nPlot = 1;
% for freqComp = [23,43,127];
%     [ftshist, binpos] = hist(spectra(freqComp,1:240));
%     subplot(1,3,nPlot); 
%     bar(binpos,ftshist,'FaceColor',[1 0 0]);hold on;
%     [ftshist, binpos] = hist(spectra(freqComp,241:end));
%     bar(binpos,ftshist,'FaceColor',[0 0 1]);hold off;
%     nPlot = nPlot+1;
% end

% Using pdf
fig = figure(11);
set(fig,'position',[365 345 875 240]);
xMax = 0;
yMax = 0;
nPlot = 1;
% Pump WDTid 32 Acc Y
% for freqComp = [24,43,127];
% Pump WDTid 32 Acc Z
for freqComp = [24,43,115];
% Pump WDTid 18 Acc Y
% for freqComp = [26,45,131];
% Pump WDTid 18 Acc Z
% for freqComp = [24,45,177];
    xValues = sort(spectra(freqComp,1:240))';
    pd = fitdist(xValues, 'Gamma');
    pdfValues = pdf(pd,xValues);
    
    if (nPlot==1)
        % Axes limits for plotting
        xMax = max(xMax,max(xValues));
        yMax = max(yMax,max(pdfValues));
    end
    
    subplot(1,3,nPlot); 
    plot(xValues,pdfValues,'r','LineWidth',2);hold on;
    % Plot the pdf corresp to the 95th percentile at this magFft
    prc95AtThisMagFft = prctileAb(freqComp);
    plot(prc95AtThisMagFft,pdf(pd,prc95AtThisMagFft),'or');
    
    xValues = sort(spectra(freqComp,241:end))';
    pd = fitdist(xValues, 'Gamma');
    pdfValues = pdf(pd,xValues);    
    
    if (nPlot==1)
        % Axes limits for plotting
        xMax = max(xMax,max(xValues));
        yMax = max(yMax,max(pdfValues));
    end
    
    plot(xValues,pdfValues,'b','LineWidth',2);
    % Plot the pdf corresp to the 95th percentile at this magFft
    prc95AtThisMagFft = prctileNo(freqComp);
    plot(prc95AtThisMagFft,pdf(pd,prc95AtThisMagFft),'ob');hold off;
    grid on;    
    xlabel('|FFT|');    
    
    title(sprintf('%d Hz',freqComp));
    axis([0 xMax+5 0 yMax+.15]);
    nPlot = nPlot+1;
end
legend('p(magFft|Ab)','p(95PrcMagFft|Ab)','p(magFft|No)','p(95PrcMgFft|No)');

if (saveFigureOption)
    % Save plot    
    plotName = fullfile(pathName,sprintf('pdfCertainFreqCompNoVsAbWdt%dAcc%s',WDTid,accelerometerSignal));
    fig = gcf;
    fig.PaperPositionMode = 'auto';    
    print(plotName,'-dpng','-r0');
    plotName = fullfile(pathName,sprintf('pdfCertainFreqCompNoVsAbWdt%dAcc%s.fig',WDTid,accelerometerSignal));
end
savefig(plotName);