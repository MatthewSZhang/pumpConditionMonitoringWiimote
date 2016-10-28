function [spectra,spd,timeStampWindow] = wiimoteRecordingsFftPerPeriod(Y,Y_hp,timeStamp,arc,locsPeriod,fileId,preprocFeatOptions)
% AM Last Modified 24/08/2016
% This function is called by the spectrawriteupAMmodifed.m script
% Does FFT per period
% Currently, the "arc" is not computed in wiimoteRecordingsPreprocess.m,
% which is required to compute speed here 
% (uncomment the corresponding lines in wiimoteRecordingsPreprocess.m to compute arc 
% and wiimoteRecordingsFftPerPeriod.m to compute speed)

% 28/10/2016 AM Modified: Parameters bundled under preprocFeatOptions
NFFT = preprocFeatOptions.NFFT;
maxClipLenAdhoc = preprocFeatOptions.maxClipLenAdhoc; % .5*preprocFeatOptions.NFFT
plotOption = preprocFeatOptions.plotOption;

fprintf('\t\tGenerating spectra...\n');

idx = 1;
locIdx = 1;  
spectra = [];
spd = [];
timeStampWindow = [];

while locIdx<length(locsPeriod)
    intstart = locsPeriod(locIdx);
    intend = locsPeriod(locIdx+1);
%     fprintf('fileId=%d, %d %d clipLen=%d clipLen<=NFFT=%d\n',fileId,intstart,intend,intend-intstart+1,intend-intstart+1<=NFFT)               

    % Skip this unusually long period (potentially due to e.g. switch
    % between pump users during the same fileIdording)
    if (intend-intstart+1>NFFT || intend-intstart+1>=maxClipLenAdhoc)
        disp('cliplength>NFFT!');        
        locIdx = locIdx+1;                    
    else
        % Farah's suggestion - Check if it's a useful pump cycle,
        %  i.e. enough periodic variation? Use original Y to check this?

        % Apply fft to the highpass filtered Y signal corresponding to one pump cycle
        Yfhat = fft(Y_hp(intstart:intend), NFFT);
        % save this up to the nyquist limit. 
        spectra(:,idx) = abs(Yfhat(1:NFFT/2));
        % Time stamp when this period starts
        timeStampWindow(idx,1) = timeStamp(intstart);
        
        if plotOption            
            figure(100);
            subplot(4,1,1);plot(Y(intstart:intend));grid on;grid minor;
            ylabel('Y');
            subplot(4,1,2);plot(Y_hp(intstart:intend));grid on;grid minor;
            ylabel('Y highpassed');
            xlabel('Freq (Hz)');
            subplot(4,1,3);plot(spectra(:,idx));grid on;grid minor;
            ylabel('|Y_fft|');            
            subplot(4,1,4);plot(10*log10(spectra(:,idx)));grid on;grid minor;
            ylabel('|Y_fft| (dB)');
            figure(100);
            pause(.05);         
        end

%         % calculate amount of movement in the window
%         spd(idx,1) = sum(abs(diff(arc(intstart:intend))));
                
        idx = idx+1;
        locIdx = locIdx+1;
    end                                
end

fprintf('\t\tSpectra generation complete\n');