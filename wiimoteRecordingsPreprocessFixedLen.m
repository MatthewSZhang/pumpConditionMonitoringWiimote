function [Y,Y_hp,arc,timeStamp] = wiimoteRecordingsPreprocessFixedLen(Y,Z,timeStamp,fileId,preprocFeatOptions)
% AM Last Modified 24/08/2016
% This function is called by the wiimoteRecordingsAnalysisScript.m script
% Does all the preprocessing
% LPF >> Find peaks >> Remove the ends of the original signal >> HPF
% Currently, the "arc" is not computed, which is required to compute speed
% (uncomment the corresponding lines to compute arc)
% Caution - Involves various adhoc numbers/parameters!

% 28/10/2016 AM Modified: Parameters bundled under preprocFeatOptions
NFFT = preprocFeatOptions.NFFT;
thresholdToFindPeaks = preprocFeatOptions.thresholdToFindPeaks; % threshold used when calling findpeaks.m
distBetnTroughsThres = preprocFeatOptions.distBetnTroughsThres; % median +/- 0.5*median
minNumOfTroughsInRecording = preprocFeatOptions.minNumOfTroughsInRecording; % In order to consider this recording
Fs = preprocFeatOptions.Fs; % sampling frequency will change depending on downsample factor
plotOption = preprocFeatOptions.plotOption;

fprintf('\t\tPreprocessing...\n');

%General gist: 
%Design of low-pass/high-pass filters:Done fairly arbitrarily, based on what looks OK

% Zero-phase filtering
% http://uk.mathworks.com/help/signal/ref/filtfilt.html
% AM Modified 27 Oct 2016: If the signal is downsampled, Fpass, Fstop should change?
Fpass = 2;% 2*Fpass(Hz)/sampling f(Hz) in normalized freq
Fstop = 4;% 2*Fstop(Hz)/sampling f(Hz) in normalized freq
Lowpass = designfilt('lowpassfir', ...
    'PassbandFrequency',2*Fpass/Fs,'StopbandFrequency',2*Fstop/Fs, ...
    'PassbandRipple',1,'StopbandAttenuation',60, ...
    'DesignMethod','equiripple');
% Lowpass = designfilt('lowpassfir', ...
%     'PassbandFrequency',2*Fpass/Fs,'StopbandFrequency',2*Fstop/Fs, ...
%     'PassbandRipple',0.5, 'StopbandAttenuation',65,...
%     'DesignMethod','kaiserwin');

% 'PassbandRipple',0.5,'StopbandAttenuation',65,'DesignMethod','kaiserwin');

% AM Modified 26 Oct 2016: Using Signal-MovingAverage(Delayed) instead
% Fpass = 4; % 2*Fpass(Hz)/sampling f(Hz) in normalized freq
% Fstop = 2; % 2*Fstop(Hz)/sampling f(Hz) in normalized freq
% Highpass = designfilt('highpassfir', ...
%     'PassbandFrequency',2*Fpass/Fs,'StopbandFrequency',2*Fstop/Fs, ...
%     'PassbandRipple',1,'StopbandAttenuation',60, ...
%     'DesignMethod','equiripple');
% Highpass = designfilt('highpassfir', ...
%     'PassbandFrequency',2*Fpass/Fs,'StopbandFrequency',2*Fstop/Fs, ...
%     'PassbandRipple',0.5,'StopbandAttenuation',65,...
%     'DesignMethod','kaiserwin');

% Visualize the filter
% fvtool(Lowpass);

% AM Modified 25 Oct 2016: Intuitively and quantitatively, shouldn’t/doesn’t affect FFT
% For testing, no need to store normalization constant from training
% Normalize (Zero mean)
% Y = Y - mean(Y);
% Z = Z - mean(Z);
arc = [];

% lowpass filter and arc calculation for gross movement
Y_lp = filtfilt(Lowpass,Y);
% Uncomment the following to compute arc
% Z_lp = filtfilt(Lowpass,Z);
% arc = atand(Y_lp./(Z_lp+Z)); % Why?

if (plotOption)
    xMax = timeStamp(end)+5;
    yMin = min(Y);
    yMax = max(Y);
end
if (plotOption)
    figure(fileId);subplot(5,1,1);plot(timeStamp,Y);ylabel('Y');grid on;
    axis([0 xMax yMin yMax]);
    title('Preprocessing');
    figure(fileId);subplot(5,1,2);plot(timeStamp,Y_lp);ylabel('Y Low-passed');grid on;
    axis([0 xMax yMin yMax]);
end

% Find troughs to find periods
% Clip ends, i.e. remove signal before/after first/last trough
x = (1:length(Y_lp))';
% Find the troughs
[~,locsTroughs] = findpeaks(-Y_lp,x,'MinPeakProminence',thresholdToFindPeaks);
if (plotOption)  
    figure(fileId);subplot(5,1,3);
    plot(timeStamp,Y_lp);hold on;
    plot(timeStamp(locsTroughs),Y_lp(locsTroughs),'or'); 
    hold off;grid on;
    axis([0 xMax yMin yMax]);
end


% There has to be at least 3 troughs!
if length(locsTroughs)<=minNumOfTroughsInRecording
    Y = [];
    Y_hp = [];
    arc = [];
    timeStamp = [];
else
    % Remove the first and last trough
    locsTroughs = locsTroughs(2:end-1);
    interval2 = (locsTroughs(1):locsTroughs(end))';
    
    % Check if the length of time-series >= NFFT
    if length(interval2)<NFFT
        Y = [];
        Y_hp = [];
        arc = [];
        timeStamp = [];
    else    
        Y = Y(interval2); % Y ends removed
        %     arc = arc(interval2); % Adjust arc
        timeStamp = timeStamp(interval2); % Adjust timeStamp        
        locsTroughs = bsxfun(@minus,locsTroughs,locsTroughs(1)-1); % Adjust locs of troughs 

        % Only consider window with "typical" observations
        % Method I: "typical" = minDist < distBetnTroughs < maxDist
        indicesToRemove = [];
        distBetnTroughs = locsTroughs(2:end)-locsTroughs(1:end-1);

        fprintf('FileId %d Pump stroke period (secs): Avg %.1f, 1SD %.1f, Med %.1f\n',...
            fileId,mean(distBetnTroughs/Fs),std(distBetnTroughs/Fs),median(distBetnTroughs/Fs));

        medDistBetnTroughs = median(distBetnTroughs);
        % Define .5*medDistBetnTroughs < typicalTroughs < 1.5*medDistBetnTroughs
        atypicalLocsTroughs = find(distBetnTroughs>(1+distBetnTroughsThres)*medDistBetnTroughs...
            | distBetnTroughs<distBetnTroughsThres*medDistBetnTroughs==1);
        for idx=1:size(atypicalLocsTroughs,1)
            indicesToRemove = cat(1,indicesToRemove,(locsTroughs(atypicalLocsTroughs(idx)):locsTroughs(atypicalLocsTroughs(idx)+1))');
        end
        Y(indicesToRemove) = [];
        % arc(indicesToRemove) = [];
        timeStamp(indicesToRemove) = [];

        % % Method II: "typical" = composed of certain no. of sinusoids
        % Y_lp = Y_lp(interval2); % Y ends removed
        % % Find the peaks
        % x = (1:length(Y_lp))';
        % [~,locsPeaks] = findpeaks(Y_lp,x,'MinPeakProminence',threshold);
        % 
        % if (plotOption)  
        %     figure(fileId);subplot(5,1,4);
        %     plot(timeStamp,Y);grid on;
        %     ylabel('Y ends removed');hold on;
        %     plot(timeStamp(locsTroughs),Y(locsTroughs),'or');
        %     plot(timeStamp(locsPeaks),Y(locsPeaks),'ok');
        %     hold off;grid on;
        %     axis([0 xMax -1 1]);
        % end    
        %
        % minCountTroughPeaks = 3;
        % maxCountTroughPeaks = 8;
        % Nwindows = floor(length(Y)/NFFT);
        % indicesToRemove = [];
        % idxEnd = 0;
        % for nWin=1:Nwindows
        %     idxStart = idxEnd+1;
        %     idxEnd = idxStart+NFFT-1;
        %     % Check if each window is "typical"
        %     countTroughsThisWindow = size(locsTroughs(locsTroughs>=idxStart & locsTroughs<=idxEnd),1);
        %     countPeaksThisWindow = size(locsPeaks(locsPeaks>=idxStart & locsPeaks<=idxEnd),1);
        %     % Following range corresponds to pumping stroke speed
        %     if (countTroughsThisWindow<minCountTroughPeaks || countTroughsThisWindow>maxCountTroughPeaks...
        %             || countPeaksThisWindow<minCountTroughPeaks && countPeaksThisWindow>maxCountTroughPeaks)
        %         fprintf('Removing indices %d : %d\n',idxStart,idxEnd);        
        %         indicesToRemove = cat(1,indicesToRemove,(idxStart:idxEnd)');        
        %     end
        % end
        % Y(indicesToRemove) = [];
        % % arc(indicesToRemove) = [];
        % timeStamp(indicesToRemove) = [];

        % Apply highpass filter to whole recording
        % AM Modified 25 Oct 2016: Using a simpler HPF with a view to
        % implement as a C-code in MPLAB 
%         % MATLAB's FIR followed by filtfilt (takes care of phase shift)
%         Y_hp = filtfilt(Highpass, Y);         
        % HPF = Signal - LPF(moving average)
        Y_hp = filterHpfFirUsingLpfMaPhaseCorrect(Y,plotOption);
        if (plotOption)
            figure(fileId);subplot(5,1,5);plot(timeStamp,Y_hp);ylabel('Y High-passed');grid on; 
            axis([0 xMax -.25 .25]);
            xlabel('Time (s)');                      
%             % Plot diff between filtfilt and HPF=LPF-MA            
%             % Moving average filter coeff   
%             figure(fileId+100);subplot(3,1,2); 
%             plot(timeStamp,Y_hp,'b');
%             hold on;            
%             plot(timeStamp,Y_hpma,'r');
%             axis([0 xMax -.4 .4]);
%             hold off;
%             legend('Yhp','Yhpma');
% 
%             figure(fileId+100);subplot(3,1,3);
%             plot(timeStamp,Y_hp-Y_hpma);
%             axis([0 xMax -.4 .4]);            
%             legend('Yhp-Yhpma');
%             suptitle('Yhp vs. Yhpma');
            drawnow;
        end        
    end
end

fprintf('\t\tPreprocessing complete\n');