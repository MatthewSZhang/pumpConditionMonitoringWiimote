function [spectra,spd,timeStampWindow] = wiimoteRecordingsSpectrogram(Y,Y_hp,timeStampThisRec,arc,fileId,preprocFeatOptions)
% AM Last Modified 24/08/2016
% This function is called by the wiimoteRecordingsAnalysisScript.m script
% Generates spectrogram using fixed length overlapping window fourier transform 
% Currently, the "arc" is not computed in wiimoteRecordingsPreprocess.m,
% which is required to compute speed here 
% (uncomment the corresponding lines in wiimoteRecordingsPreprocessFixedLen.m to compute arc 
% and wiimoteRecordingsSpectrogram.m to compute speed)

% 28/10/2016 AM Modified: Parameters bundled under preprocFeatOptions
NFFT = preprocFeatOptions.NFFT;
overlapFrac = preprocFeatOptions.overlapFrac; % .5*preprocFeatOptions.NFFT

fprintf('\t\tGenerating spectra...\n');

if length(Y_hp)<NFFT
    spd = [];
    spectra = [];    
    timeStampWindow = [];
else
    spd = [];
    % Spectrogram of accelerometry data
    nsc = NFFT;                      % lenth of each segment
    nov = floor(overlapFrac*NFFT);   % No of samples overlap
    nff = NFFT;                      % length of FFT
    spectra = spectrogram(Y_hp,nsc,nov,nff,'yaxis');
    
    % Find approximate corresponding time stamps
    Nwindows = size(spectra,2);    
    timeStampWindow = linspace(timeStampThisRec(NFFT/2),timeStampThisRec(end-NFFT/2),Nwindows)';    
end

fprintf('\t\tSpectra generation complete\n');