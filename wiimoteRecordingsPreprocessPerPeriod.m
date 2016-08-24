function [Y,Y_hp,arc,locsTroughs,timeStamp] = wiimoteRecordingsPreprocessPerPeriod(Y,Z,timeStamp,fileId,plotOption)
% AM Last Modified 24/08/2016
% This function is called by the wiimoteRecordingsAnalysisScript.m script
% Does all the preprocessing
% LPF >> Find peaks >> Remove the ends of the original signal >> HPF
% Currently, the "arc" is not computed, which is required to compute speed
% (uncomment the corresponding lines to compute arc)
% Caution - Involves various adhoc numbers/parameters!

fprintf('\t\tPreprocessing...\n');

thresholdToFindPeaks = .25; % threshold used when calling findpeaks.m
minNumOfTroughsInRecording = 3; % In order to consider this recording

%General gist: 
%Design of high pass filte:Done fairly arbitrarily, based on what looks OK
Fs = 96;  % Sampling Frequency

% Zero-phase filtering
% http://uk.mathworks.com/help/signal/ref/filtfilt.html
Fpass = 2; % 2*Fpass(Hz)/sampling f(Hz) in normalized freq
Fstop = 4; % 2*Fstop(Hz)/sampling f(Hz) in normalized freq
Lowpass = designfilt('lowpassfir', ...
    'PassbandFrequency',2*Fpass/Fs,'StopbandFrequency',2*Fstop/Fs, ...
    'PassbandRipple',1,'StopbandAttenuation',60, ...
    'DesignMethod','equiripple');
% Lowpass = designfilt('lowpassfir', ...
%     'PassbandFrequency',2*Fpass/Fs,'StopbandFrequency',2*Fstop/Fs, ...
%     'PassbandRipple',0.5, 'StopbandAttenuation',65,...
%     'DesignMethod','kaiserwin');

% 'PassbandRipple',0.5,'StopbandAttenuation',65,'DesignMethod','kaiserwin');

Fpass = 4; % 2*Fpass(Hz)/sampling f(Hz) in normalized freq
Fstop = 2; % 2*Fstop(Hz)/sampling f(Hz) in normalized freq
Highpass = designfilt('highpassfir', ...
    'PassbandFrequency',2*Fpass/Fs,'StopbandFrequency',2*Fstop/Fs, ...
    'PassbandRipple',1,'StopbandAttenuation',60, ...
    'DesignMethod','equiripple');
% Highpass = designfilt('highpassfir', ...
%     'PassbandFrequency',2*Fpass/Fs,'StopbandFrequency',2*Fstop/Fs, ...
%     'PassbandRipple',0.5,'StopbandAttenuation',65,...
%     'DesignMethod','kaiserwin');

% Visualize the filter
% fvtool(Lowpass);

Y = Y - mean(Y);
Z = Z - mean(Z);
arc = [];

if length(Y) > 396
    % lowpass filter and arc calculation for gross movement
    Y_lp = filtfilt(Lowpass,Y);
    Z_lp = filtfilt(Lowpass,Z);
%     arc = atand(Y_lp./(Z_lp+Z)); % Why?

    if (plotOption)
        xMax = timeStamp(end)+5;
        figure(fileId);subplot(5,1,1);plot(timeStamp,Y);ylabel('Y');grid on;
        axis([0 xMax -1 1]);
        title('Preprocessing');
        figure(fileId);subplot(5,1,2);plot(timeStamp,Y_lp);ylabel('Y Low-passed');grid on;
        axis([0 xMax -1 1]);
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
        axis([0 xMax -1 1]);
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
        Y = Y(interval2); % Y ends removed
    %     arc = arc(interval2); % Adjust arc
        timeStamp = timeStamp(interval2); % Adjust time        
        locsTroughs = bsxfun(@minus,locsTroughs,locsTroughs(1)-1); % Adjust locsPeriod of troughs 
        if (plotOption)  
            figure(fileId);subplot(5,1,4);plot(timeStamp,Y);grid on;
            ylabel('Y ends removed');
            axis([0 xMax -1 1]);
        end        

        % Apply highpass filter to whole recording
        Y_hp = filtfilt(Highpass, Y); 
        if (plotOption)
            figure(fileId);subplot(5,1,5);plot(timeStamp,Y_hp);ylabel('Y High-passed');grid on; 
            axis([0 xMax -.25 .25]);
            xlabel('Time (s)');
        end 
    end
end

fprintf('\t\tPreprocessing complete\n');