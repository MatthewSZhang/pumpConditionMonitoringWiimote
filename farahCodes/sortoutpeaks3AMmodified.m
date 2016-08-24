function [locs,pks] = sortoutpeaks3AMmodified(Y_s, threshold)
% AM Modified 30/06/2016
% This function is called by the wiimoteRecordingsPreprocess.m script

x = (1:length(Y_s))';

% Find the troughs
[pks,locs] = findpeaks(-Y_s,x,'MinPeakProminence',threshold);