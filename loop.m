%This script will loop through all the csv files in the current folder and
%run ReadNormSSCDRCData.m then hillFitv2.m

%be sure to modify endRow and the filename to fit your needs.  Every
%replicate must be in a different CSV file in order for the scripts to work
%properly.

concentrations = [10.000000
2.000000
0.400000
0.080000
0.016000
0.003200
0.000640
0.000128
5.000000
1.000000
0.200000
0.040000
0.008000
0.001600
0.000320
0.000064];

bigstructNormed = struct;
%this is going to be a nested struct containing ALL the data.  First level
%is different replicate days.  Second level is separate plates on each day.
list=dir('/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/DRCs on Resistant Cells/Resistant Cells Round 3/IQueCSVOutput/*.csv');
for i = 1:length(list);
    filename = sprintf('%s','/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/DRCs on Resistant Cells/Resistant Cells Round 3/IQueCSVOutput/',list(i).name);
    endRow = 5181;
    bigstructNormed.((sprintf('Rep%d',i))) = ReadNormSSCDRCData(filename,endRow);
end

%create new struct, name fields by looking at fields of structs already
%have, stripping date, and checking if there's already a field with that
%name.  Put data from that field horizcat into the new field.

concatinatedStruct = struct;
numReps = length(list);
for sepStruct = 1:numReps
    fields = fieldnames(bigstructNormed.((sprintf('Rep%d',sepStruct))));
    for field = 1:length(fields) 
        fieldName = char(fields(field,:));
        newFieldName = fieldName(13:end);
        if isfield (concatinatedStruct,newFieldName)
            %horizcat data into whats already there
            together = horzcat(concatinatedStruct.(newFieldName),bigstructNormed.((sprintf('Rep%d',sepStruct))).(char(fieldName)));
            concatinatedStruct.(newFieldName) = together;
        else 
            %make a new field with that name and put the data in there
            concatinatedStruct.(newFieldName) = bigstructNormed.((sprintf('Rep%d',sepStruct))).(char(fieldName));
        end
    end
end
    
[fittedStruct,individualFittedStruct] = hillFitv2(concatinatedStruct,concentrations,numReps);

%Export all IC50s into a matrix
cellLines = fieldnames(fittedStruct);
drugNames = {'Doxorubicin'
    'Vincristine'
    'YM155'
    'Simvastatin'
    'Fluvistatin'
    'Camptothecan'
    '5-FU'
    'Actinomycin D'
    'SAHA'
    'VER-50589'
    'Etopiside'
    'Cisplatin'
    'Colchicine'
    'Gemcytabine'
    'Paclitaxel'
    'Puromycin'
    'Scriptaid'
    '17-AAG'};
IC50s = zeros(length(cellLines),length(drugNames));
for field = 1:length(cellLines)
    IC50s(field,:) = fittedStruct.(char(cellLines(field)))(3,3:21);
end

%Calculate log2 fold change over eu-myc.  First find eumyc!
euMycIndex = find(~cellfun(@isempty,strfind(cellLines,'Myc')));
log2FoldChanges = zeros(size(IC50s));
for drug = 1:length(drugNames)
    log2FoldChanges(:,drug) = log2(IC50s(:,drug)/IC50s(euMycIndex,drug));
end

%Remove all the eu-myc ones
log2FoldChangesnoEu = [log2FoldChanges(1:euMycIndex-1,:);log2FoldChanges(euMycIndex+1:end,:)];
cellLinesnoEu = vertcat(cellLines(1:euMycIndex-1,:),cellLines(euMycIndex+1:end,:));

%Get values for heatmap scale
log2ScaleMax = max(max(log2FoldChangesnoEu));
log2ScaleMin = min(min(log2FoldChangesnoEu));
log2Scale = max(log2ScaleMax,abs(log2ScaleMin));

%Make heatmap!
foldChangeHeatmap = clustergram(log2FoldChangesnoEu','RowLabels',drugNames,'ColumnLabels',cellLinesnoEu,'DisplayRange',log2Scale,'Symmetric','true','Colormap',redbluecmap,'ColumnLabelsRotate',0);
plot(foldChangeHeatmap);
colormap(redbluecmap(256));
imagesc(log2FoldChangesnoEu',[-log2Scale log2Scale]);
colorbar;