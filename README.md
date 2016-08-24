# pumpConditionMonitoringWiimote
On-board analysis code (preprocessing, feature generation, and classification)

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

