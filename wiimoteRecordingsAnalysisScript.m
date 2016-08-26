clear
close all

% AM Last Modified 24/08/2016
%
% Notes (PLEAE READ THESE NOTES!)
% Comment/uncomment the start of for loop corresponding to the WDT ID
% Set the plotOption = true or false, depending on whether you want plots
% Regarding labels, currently only available for WDT ID 18 and 32 (please
% uncomment the corresponding lines below to define/plot these labels)

% Mostly based on Farah's codes @Z:\netshares\bspprojects6\SmartWater\Scripts\June16
% Differences
% Filters designed using MATLAB's designfilt function
% wiimoteRecordingsPreprocess***.m function does all the preprocessing
% wiimoteRecordings***.m function does FFT per period
% fileIdLookupTableFun.m function looks up fileIds corresponding to WDTid from the XLS spreadsheet

% Choose which accelerometer signal to use (Y or Z)

% Choose which preprocessing technique to use:
% (I) wiimoteRecordingsPreprocessPerPeriod.m - Preprocesses signal to do FFT per period
% (II) wiimoteRecordingsPreprocessFixedLen.m - Preprocesses signal to do fixed length FFT with overlapping window

% Choose which method to use to generate spectra:
% (I) wiimoteRecordingsFftPerPeriod - FFT per period
% (II) wiimoteRecordingsFftFixedLen.m - FFT fixed length with overlapping windows
% (III) wiimoteRecordingsSpectrogram.m - Spectrogram with overlapping windows

% Feature generation (adhoc): 
% (I) Uniformly samples across freq range
% (II) Handpicked looking at spectra plot
% Please uncomment the last part of the script to generate features and
% create/save a dataset that will be used by ndtoolboxWorkspace.m script to
% perform classification

% For each window in the clip, size of windows and their overlap = 
% Trade off between frequency/time resolution and smoothness
NFFT = 512; % At Fs = 96 Hz and Median pump stroke period = 1.2s, 512 samples = 4.5 pump strokes
maxPeriodLenAdhoc = 150; % When spectra generated using wiimoteRecordingsPreprocessPerPeriod.m
overlapFrac = .5; % .5*NFFT, When spectra generated using wiimoteRecordingsFftFixedLen.m

% Path to data
rawDataPath = 'C:\Users\engs1602\research\data\Data Feb-Mar 2016\Raw Data - wiimote';

% Choose which accelerometer signal to use 
% accelerometerSignal = 'Y';
accelerometerSignal = 'Z'; 

% Choose which preprocessing technique to use:
% (I) wiimoteRecordingsPreprocessPerPeriod.m - Preprocesses signal to do FFT per period
% (II) wiimoteRecordingsPreprocessFixedLen.m - Preprocesses signal to do fixed length FFT with overlapping window
% preprocessMethod = 'PerPeriod';
preprocessMethod = 'FixedLen';

% Choose which method to use to generate spectra:
% (I) wiimoteRecordingsFftPerPeriod - FFT per period
% (II) wiimoteRecordingsFftFixedLen.m - FFT fixed length with overlapping windows
% (III) wiimoteRecordingsSpectrogram.m - Spectrogram with overlapping windows
% spectraGenMethod = 'FftPerPeriod';
spectraGenMethod = 'FftFixedLen';
% spectraGenMethod = 'Spectrogram';

% Save figure options (Will need to modify path below as required)
saveFigureOption = false;

% Begin code
% Initialize following to help with plotting later
WDTidVec = [];
fileIdVec = [];
idx = 1;
spectra = [];
timeStampWindow = [];

% The following are condition labels based on field notes. Only to be used
% as a guideline. Not groundtruth!
conditionLabels = {'WorkingExcellent','NoisyButWorking','DryBorehole',...
    'RisingMainLeak(Rare)','BrokenSeal(VeryCommon)','Unsure',...
    'WaterLeakingFromPump','WornBushBearing','StiffHandle','WornSeal'};
    
% Condition 1: Working excellent 
% condition = 1; WDTids = [61, 133, 243, 171, 185, 117, 27];
% % Condition 2: Noisy but working
% condition = 2; WDTids = 18%[273, 18, 140]; % one more private (no WDT)
% % Condition 3: Dry borehole
% condition = 3; WDTids = [145,21];
% % Condition 4: Rising main leak (Rare)
% condition = 4; WDTids = [18];
% % Condition 5: Broken seal (Most common)
condition = 5; WDTids = [32];  
% % Condition 6: Unsure
% condition = 6; WDTids = [129,132]; % one more India Mark II (no WDT)
% % Condition 7: Water leaking from pump
% condition = 7; WDTids = [112, 262, 126, 204];
% % Condition 8: Worn bush bearing
% condition = 8; WDTids = [196, 125];
% % Condition 9: Stiff handle
% condition = 9; WDTids = [134, 97, 1]; % one more Ukunda (not in pilot, no WDT)
% % Condition 10: Worn seal
% condition = 10; WDTids = [192];   
for WDTid = WDTids
    fprintf('Analyzing data from WDT ID %d ...\n',WDTid);
    % look up fileIds corresponding to WDTid from the XLS spreadsheet
    fileIdsThisWDT = fileIdLookupTableFun(WDTid)';
    disp(fileIdsThisWDT);
    if ~any(isnan(fileIdsThisWDT))
        Yall = [];
        Zall = [];
        fileIdIndices= [];
        for fileId=fileIdsThisWDT
            fprintf('\t File ID %d\n',fileId);
            file = fopen(fullfile(rawDataPath,sprintf('%d.txt',fileId)),'rt');
            raw = textscan(file, '%f:%f:%f:%f \t %f \t %f \t %f \t %f \t %f \t');
            fclose(file);
            timeStamp = raw{4}/1000 + raw{3} + 60*raw{2};
            Y = raw{6}; 
            Z = raw{7};          
       
            % Preprocess:
            % LPF >> Find peaks >> Remove the ends of the original signal >> HPF
            % Currently, not using the output "arc" to compute speed.    
            plotOption = false;
            % Y accelerometer signal
            if strcmp(accelerometerSignal,'Y')                
                if strcmp(preprocessMethod,'PerPeriod')
                    [Y,Y_hp,arc,locsPeriod,timeStampThisRec] = wiimoteRecordingsPreprocessPerPeriod(Y,Z,timeStamp,fileId,plotOption);
                elseif strcmp(preprocessMethod,'FixedLen')
                    [Y,Y_hp,arc,timeStampThisRec] = wiimoteRecordingsPreprocessFixedLen(Y,Z,timeStamp,fileId,NFFT,plotOption);
                end
            % Z accelerometer signal
            elseif strcmp(accelerometerSignal,'Z')    
                if strcmp(preprocessMethod,'PerPeriod')
                    [Z,Z_hp,arc,locsPeriod,timeStampThisRec] = wiimoteRecordingsPreprocessPerPeriod(Z,Y,timeStamp,fileId,plotOption);
                elseif strcmp(preprocessMethod,'FixedLen')
                    [Z,Z_hp,arc,timeStampThisRec] = wiimoteRecordingsPreprocessFixedLen(Z,Y,timeStamp,fileId,NFFT,plotOption);
                end
            end
%             keyboard;

            % Generate Spectra:
            % Currently, not using the input "arc", which is required to compute speed.
            % Consequently, not getting the output "spd"
            % "timeStampWindow" = time stamp may not correspond to true time stamp due to preprocessing (Need to verify)
            plotOption = false;
            % Y accelerometer signal
            if strcmp(accelerometerSignal,'Y')
                if strcmp(spectraGenMethod ,'FftPerPeriod')
                    [spectraThisRec,spdThisRec,timeStampWindowThisRec] = wiimoteRecordingsFftPerPeriod(Y,Y_hp,timeStampThisRec,arc,locsPeriod,NFFT,maxPeriodLenAdhoc,fileId,plotOption);        
                elseif strcmp(spectraGenMethod ,'FftFixedLen')
                    [spectraThisRec,spdThisRec,timeStampWindowThisRec] = wiimoteRecordingsFftFixedLen(Y,Y_hp,timeStampThisRec,arc,NFFT,overlapFrac,fileId,plotOption);        
                elseif strcmp(spectraGenMethod ,'Spectrogram')
                    [spectraThisRec,spdThisRec,timeStampWindowThisRec] = wiimoteRecordingsSpectrogram(Y,Y_hp,timeStampThisRec,arc,NFFT,overlapFrac,fileId,plotOption);        
                end
            % Z accelerometer signal
            elseif strcmp(accelerometerSignal,'Z')        
                if strcmp(spectraGenMethod ,'FftPerPeriod')
                    [spectraThisRec,spdThisRec,timeStampWindowThisRec] = wiimoteRecordingsFftPerPeriod(Z,Z_hp,timeStampThisRec,arc,locsPeriod,NFFT,maxPeriodLenAdhoc,fileId,plotOption);
                elseif strcmp(spectraGenMethod ,'FftFixedLen')
                    [spectraThisRec,spdThisRec,timeStampWindowThisRec] = wiimoteRecordingsFftFixedLen(Z,Z_hp,timeStampThisRec,arc,NFFT,overlapFrac,fileId,plotOption); 
                elseif strcmp(spectraGenMethod ,'Spectrogram')
                    [spectraThisRec,spdThisRec,timeStampWindowThisRec] = wiimoteRecordingsSpectrogram(Z,Z_hp,timeStampThisRec,arc,NFFT,overlapFrac,fileId,plotOption);  
                end
            end
            lenThisFile = size(spectraThisRec,2);

            spectra = cat(2,spectra,spectraThisRec);
            fileIdVec = cat(1,fileIdVec,fileId*ones(lenThisFile,1));
            WDTidVec = cat(1,WDTidVec,WDTid*ones(lenThisFile,1));
            timeStampWindow = cat(1,timeStampWindow,timeStampWindowThisRec);
        end        
    end
end

%%
% Once the spectra is generated above, plot the spectra, etc.
plotOption = true;
if plotOption
    % Labels - Based on prefix/postfix
    % Prefix/postfix labels are not entirely indicative of pump working/not
    % working. Pump working/not-working is also sorta dependent case basis.
    % Better to confirm with Heloise/Farah
    isnor = ones(size(spectra,2),1);
    % WDT 18    
    % isnor(fileIdVec<74) = 0;
    % WDT 32
    isnor(fileIdVec<190) = 0;

    freqStart = 1; 
    freqEnd = 256; 
    
%     % Fixed length overlapping window FFT (STFT)
%     % Plot spectra, fileId indices, and labels
%     spectraClipped = spectra(freqStart:end,:);    
%     figure(1);
%     subplot(4,1,1:2);imagesc(spectraClipped);
%     caxis(prctile(spectraClipped(:), [5,95]))
%     % caxis([0 2]);%colorbar;
%     set(gca,'YDir','normal'); ylabel('Frequency (Hz)'); 
%     title(sprintf('|FFT| (%s)',conditionLabels{condition}));
    
    % Fixed length overlapping window FFT (STFT)
    % Plot spectra (dB), fileId indices, and labels
    spectraDB = 10*log10(spectra);
    spectraClippedDB = spectraDB(freqStart:freqEnd,:);
    % figure;hist(spectraDB1,500);
    figure(2);
    subplot(5,1,1:2);imagesc(spectraClippedDB); 
%     caxis(prctile(spectraClippedDB(:), [0,100]))
    caxis([-1 8]);%colorbar;
    set(gca,'YDir','normal'); ylabel('Frequency'); 
    title(sprintf('|FFT| dB (%s)',conditionLabels{condition}));
    
%     % Sectrogram
%     magSpectra = abs(spectra);
%     figure(2);
%     subplot(4,1,1:2);imagesc(abs(spectra));
%     caxis(prctile(magSpectra(:), [5,95]))
%     set(gca,'YDir','normal'); ylabel('Frequency (Hz)'); 
%     title(sprintf('|Spectrogram| (%s)',conditionLabels{condition})); 
    
    subplot(5,1,3);
    yyaxis left;
    plot(WDTidVec,'LineWidth',1.5); 
    ylabel('WDTid');    
    axis([0 size(spectra,2) min(WDTidVec)-1 max(WDTidVec)+1]); 
    yyaxis right;
    plot(fileIdVec,'LineWidth',1.5);grid on;grid minor;
    ylabel('FileId');   
    axis([0 size(spectra,2) min(fileIdVec)-1 max(fileIdVec)+1]);
%     minYaxis = min(cat(1,fileIdVec,WDTids'));
%     maxYaxis = max(cat(1,fileIdVec,WDTids'));
%     axis([0 size(spectra,2) minYaxis-1 maxYaxis+1]); 
%     legend('WDTid','FileId');
    subplot(5,1,4);plot(timeStampWindow,'LineWidth',1.5);ylabel('Time (s)');grid minor;axis([0 size(spectra,2) 0 max(timeStampWindow)]); 
    subplot(5,1,5);plot(isnor,'LineWidth',1.5);grid on;axis([0 size(spectra,2) -0.1 1.1]); 
%     ylabel('Prefix/postfix label');
%     xlabel('Period count across all recordings for Pump 32');
    ylabel('Prefix/postfix label');
    xlabel('Window count across all recordings');    
    
%     if (saveFigureOption)
%         % Save plot
%         pathName = fullfile('C:\Users\engs1602\research\meetings\smallGroup\20160825ManandharAnalysisCodeUpdate\plots',sprintf('caseStudyPumpWDTid%dAcc%szmuv',WDTid,accelerometerSignal));
%         plotName = fullfile(pathName,sprintf('spectraFftFixedLenWdtid%dAcc%s',WDTid,accelerometerSignal));
%         fig = gcf;
%         fig.PaperPositionMode = 'auto';    
%         print(plotName,'-dpng','-r0');
%         plotName = fullfile(pathName,sprintf('spectraFftFixedLenWdtid%dAcc%s.fig',WDTid,accelerometerSignal));
%         savefig(plotName);
%     end
end

%% 
% Once the spectra is generated above, choose freq components as features
% Adhoc: 
%     (1) Uniformly samples across freq range
%     (2) Handpicked looking at spectra plot
%% WDT 18
% Uniformly down sample along freq
% data.x = spectra(11:5:end,:)';
% Handpick interesting freq
% data.x = spectra(90:95,:)';
% data.y = isnor;
% data.name = {'WDT18RecPartial_end'};
% save('C:\Users\engs1602\research\data\pumpWiimoteND\WDT18RecPartial11Freq90_95','data');
% save('C:\Users\engs1602\research\data\pumpWiimoteND\WDT18RecPartial1111_10_end','data');
% save('C:\Users\engs1602\research\data\pumpWiimoteND\WDT18RecPartial1194_95','data');
%% WDT 32
% Uniformly down sample along freq, FFT per period
% data.x = spectraDB(56:2:64,:)';
% Hand-pick based on the median/mean plot comparison
data.x = spectraDB([24,43,115],:)';
% Hand-pick (lower-freq), FFT per period
% data.x = spectraDB(12:3:24,:)';
% Hand-pick (higher-freq), Fixed-length (512) FFT 
% data.x = spectraDB(111:4:130,:)';
data.y = isnor;
data.name = {sprintf('WDT32Rec%s_Freq24_43_115dB',accelerometerSignal)};
save(fullfile('C:\Users\engs1602\research\data\pumpWiimoteND\expFftHandPickZmeanVsZmuv\',sprintf('WDT32Rec%s_Freq24_43_115dBzmean',accelerometerSignal)),'data');
% 
% figure(3);
% subplot(3,1,1:2);imagesc(data.x');caxis(prctile(spectraClippedDB(:), [5,95]))
% subplot(3,1,3);plot(data.y);axis([0 size(data.x,1) 0 1]);
