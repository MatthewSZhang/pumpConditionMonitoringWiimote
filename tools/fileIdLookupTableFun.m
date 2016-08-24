function fileIdsThisWDT = fileIdLookupTableFun(WDTidQuery)

fieldNotesPath = 'C:\Users\engs1602\research\data\Data Feb-Mar 2016\16.03.14 Field data notes.xlsx';
sheet = 3;
xlRange = 'A2:A314';
fileId = xlsread(fieldNotesPath,sheet,xlRange);
xlRange = 'F2:F314';
WDTid = xlsread(fieldNotesPath,sheet,xlRange);

fileIdsThisWDT = fileId(WDTid==WDTidQuery)