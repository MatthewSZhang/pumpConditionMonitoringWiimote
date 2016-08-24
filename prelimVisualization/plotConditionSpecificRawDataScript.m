clear
close all

% AM Last Modified 09/08/2016
%
% Notes (PLEAE READ THESE NOTES!)
% Comment/uncomment the start of for loop corresponding to the WDTids
% Set the plotOption = true or false, depending on whether you want plots
            
plotOption = true;

rawDataPath = 'C:\Users\engs1602\research\data\Data Feb-Mar 2016\Raw Data - wiimote';
% Choose which accelerometer signal to use (Y or Z)
accelerometerSignal = 'Y'; % 'Z'

% Condition 1: Working excellent 
% WDTids = [61, 133, 243, 171, 185, 117, 27];
% Condition 2: Noisy but working
% WDTids = [273];
% Condition 3: Dry borehole
% WDTids = [21];
% Condition 4: Rising main leak (Rare)
% WDTids = [18];
% Condition 5: Broken seal (Most common)
% WDTids = [32];  
% Condition 6: Unsure
% Condition 7: Water leaking from pump
WDTids = [32, 112, 262, 126, 204];
% Condition 8: Worn bush bearing
% WDTids = [196];
% Condition 9: Stiff handle
% WDTids = [134, 97, Ukunda (not in pilot, no WDT), 1]
% Condition 10: Worn seal
% WDTids = [192];   
for WDTid = WDTids
    fileIdsThisWDT = fileIdLookupTableFun(WDTid)';
    if ~any(isnan(fileIdsThisWDT))
        Yall = [];
        Zall = [];
        fileIdIndices= [];
        for fileId=fileIdsThisWDT
            file = fopen(fullfile(rawDataPath,sprintf('%d.txt',fileId)),'rt');
            raw = textscan(file, '%f:%f:%f:%f \t %f \t %f \t %f \t %f \t %f \t');
            fclose(file);
%             time = raw{4}/1000 + raw{3} + 60*raw{2};
            Y = raw{6}; 
            Z = raw{7};          
            
            % Normalize (Zero mean)
            Y = Y - mean(Y);
            Z = Z - mean(Z);

            Yall = cat(1,Yall,Y);
            Zall = cat(1,Zall,Z);
            lenThisFile = size(Y,1);
            fileIdIndices = cat(1,fileIdIndices,fileId*ones(lenThisFile,1));
        end
        if plotOption
            lenAllFilesThisWdtid = size(fileIdIndices,1);
            figure(WDTid);            
            subplot(3,1,1);plot(Yall);hold on;                        
            axis([0 lenAllFilesThisWdtid -2 2]);
            subplot(3,1,2);plot(Zall);
            axis([0 lenAllFilesThisWdtid -2 2]);
            subplot(3,1,3);plot(fileIdIndices);
            axis([0 lenAllFilesThisWdtid min(fileIdIndices) max(fileIdIndices)]);
        end
    end    
end