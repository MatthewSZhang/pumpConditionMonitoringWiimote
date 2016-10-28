function [spectra,spd,timeStampWindow] = wiimoteRecordingsFftFixedLen(Y,Y_hp,timeStampThisRec,arc,fileId,preprocFeatOptions)
% AM Last Modified 24/08/2016
% This function is called by the wiimoteRecordingsAnalysisScript.m script
% Performs fixed length FFT with overlapping window (STFT)
% Currently, the "arc" is not computed in wiimoteRecordingsPreprocess.m,
% which is required to compute speed here 
% (uncomment the corresponding lines in wiimoteRecordingsPreprocessFixedLen.m to compute arc 
% and wiimoteRecordingsFftFixedLen.m to compute speed)

% 28/10/2016 AM Modified: Parameters bundled under preprocFeatOptions
NFFT = preprocFeatOptions.NFFT;
overlapFrac = preprocFeatOptions.overlapFrac; % .5*preprocFeatOptions.NFFT
plotOption = preprocFeatOptions.plotOption;

fprintf('\t\tGenerating spectra...\n');

% AM Modified 25 Oct 2016: Simulate input/output of MPLAB C-code FFT
% implementation, i.e. 
% Input to FFT = round(scaled by 1000)?
% Output of FFT = round(mag(fft))

if length(Y_hp)<NFFT
    spd = [];
    spectra = [];    
    timeStampWindow = [];
else
    % FFT of accelerometry data overlapping windows
    spd = [];
    Nwindows = length((1:1-overlapFrac:floor(length(Y_hp)/NFFT))');
    % Find approximate corresponding time stamps    
    timeStampWindow = linspace(timeStampThisRec(NFFT/2),timeStampThisRec(end-NFFT/2),Nwindows)';        
    spectra = nan(NFFT/2,Nwindows);
        
    idx = 1;
    for nWin=1:Nwindows
        idxStart = (idx-1)*overlapFrac*NFFT+1;
        idxEnd = idxStart+NFFT-1;
    %     fprintf('%d %d %d\n',idx,idxStart,idxEnd);           
        
        % Apply FFT per window to the highpassed Y signal         
        Yfhat = fft(Y_hp(idxStart:idxEnd), NFFT);
        
        % save this up to the nyquist limit. 
        spectra(:,nWin) = abs(Yfhat(1:NFFT/2));
    
        if plotOption            
            figure(100);
            subplot(4,1,1);plot(Y(idxStart:idxEnd));grid on;grid minor;
            ylabel('Y');
            subplot(4,1,2);plot(Y_hp(idxStart:idxEnd));grid on;grid minor;
            ylabel('Y highpassed');
            xlabel('Freq (Hz)');
            subplot(4,1,3);plot(spectra(:,nWin));grid on;grid minor;
            ylabel('|Y_fft|');            
            subplot(4,1,4);plot(10*log10(spectra(:,nWin)));grid on;grid minor;
            ylabel('|Y_fft| (dB)');
            figure(100);
            pause(.05);         
        end
    
    %         % calculate amount of movement in the window
    %         spd(nWin,1) = sum(abs(diff(arc(intstart:intend))));      
        idx = idx+1;
    end
end

fprintf('\t\tSpectra generation complete\n');