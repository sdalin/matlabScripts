%This script will loop through all the txt files in the current folder and
%run ReadNormSSCDRCData.m then hillFitv2.m

concentrations = [10.000000
2.000000
0.400000
0.080000
0.016000
0.003200
0.000640
0
5.000000
1.000000
0.200000
0.040000
0.008000
0.001600
0.000320
0];
list=dir('/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/SSC Heterogeneity DRCs/All CSV Files/*.txt');
for i = 1:length(list);
    filename = sprintf('%s','/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/SSC Heterogeneity DRCs/All CSV Files/160511_SSC_DRC_Fixed.csv',list(i).name);
    if i == 1
        endRow = 1211;
    else
        endRow = 1933;
    end
    [bigstructNormed] = ReadNormSSCDRCData(filename,endRow);
    [fittedStruct] = hillFitv2(bigstructNormed,concentrations);
end
