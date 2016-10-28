function [X,timeStamp] = wiimoteRecordingsDownsampleByMean(X,timeStamp,preprocFeatOptions)
% Downsample by using mean of consecutive samples

downsampleFactor = preprocFeatOptions.downsampleFactor;
plotOption = preprocFeatOptions.plotOption;

remainder = rem(size(X,1),downsampleFactor);
if remainder~=0
    X = X(remainder+1:end);
    timeStamp = timeStamp(remainder+1:end);
end

if plotOption
    figure(1);
    plot(timeStamp,X,'b');
end

patternRepeatTimes = downsampleFactor;
patternTotalLength = size(X,1)/patternRepeatTimes;                        
patternIndicesTemp = repmat(1:patternTotalLength,[patternRepeatTimes 1]);
patternIndices = patternIndicesTemp(:);
clear patternRepeatTimes patternTotalLength patternIndicesTemp;
X = accumarray(patternIndices,X,[],@(x) mean(x));            
timeStamp = accumarray(patternIndices,timeStamp,[],@(x) mean(x));

if plotOption
    figure(1);
    hold on;plot(timeStamp,X,'r');hold off;
    legend('Raw','Downsampled');
    title('Downsampled by taking mean of consecutive samples');
end
            
            