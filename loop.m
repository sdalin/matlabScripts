%This script will loop through all the csv files in the current folder and
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

bigstructNormed = struct;
%this is going to be a nested struct containing ALL the data.  First level
%is different replicate days.  Second level is separate plates on each day.
list=dir('/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/SSC Heterogeneity DRCs/All CSV Files/*.csv');
for i = 1:length(list);
    filename = sprintf('%s','/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/SSC Heterogeneity DRCs/All CSV Files/',list(i).name);
    endRow = 1933;
    bigstructNormed.((sprintf('Rep%d',i))) = ReadNormSSCDRCData(filename,endRow);
end

%create new struct, name fields by looking at fields of structs already
%have, stripping date, and checking if there's already a field with that
%name.  Put data from that field vertcat into the new field.

concatinatedStruct = struct;
for sepStruct = 1:length(list)
    fields = fieldnames(bigstructNormed.((sprintf('Rep%d',sepStruct))));
    for field = 1:length(fields) 
        fieldName = char(fields(field,:));
        newFieldName = fieldName(13:end);
        if isfield (concatinatedStruct,newFieldName)
            %vertcat data into whats already there
            together = vertcat(concatinatedStruct.(newFieldName),bigstructNormed.((sprintf('Rep%d',sepStruct))).(char(fieldName)));
            concatinatedStruct.(newFieldName) = together;
        else 
            %make a new field with that name and put the data in there
            concatinatedStruct.(newFieldName) = bigstructNormed.((sprintf('Rep%d',sepStruct))).(char(fieldName));
        end
    end
end
    

[fittedStruct] = hillFitv2(bigstructNormed,concentrations);