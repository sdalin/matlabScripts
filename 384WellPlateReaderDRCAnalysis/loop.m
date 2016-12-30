%This script will loop through all the csv files in the current folder and
%run ReadNormPlateReaderDRCData.m then hillFitv2.m

%Data format must be from Tecan Plate Reader in order to work.  Current
%format is designed to read dose responses with different drugs in
%different columns, no cells in column 23, and no drugs in columns 1, 2,
%and 24.  Dose responses go from [high] in row A down, in four steps.  

%Data collected is absorbance of rezasurin salt.  So normalization first
%subtracts median of no cells wells, then divides by median of no drug
%wells.

%be sure to modify the filename to fit your needs.

load('/Users/sdalin/Dropbox (MIT)/Biology PhD/Matlab Scripts/384wellPlateReaderDRCAnalysis/dataAfterFit');
bigstructNormed = struct;
%this is going to be a nested struct containing ALL the data.  First level
%is different days of data collection.  Second level is separate plates on each day and one plate with 'info' ie, drug names and concentrations.
folderName = '/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/DRCs on Resistant Cells/Resistant Cells Round 3/DRCs 161115/Data/';
list=dir(sprintf('%s*.xlsx',folderName));
for i = 1:length(list);
    filename = sprintf('%s%s',folderName,list(i).name);
    bigstructNormed.(sprintf('date%s',filename(end-14:end-9))) = ReadNormPlateReaderDRCData(filename);
end

%Get fits of all raw data in this folder  
if ~isdir(sprintf('%s../matlabOutput',folderName))
    mkdir(sprintf('%s../matlabOutput',folderName));
end
folderName = sprintf('%s..',folderName);

%Calculate fits of individual data, all replicates, and put graphs pf individual fits in the
%folder above.
[dataAfterFit] = hillFitv2(bigstructNormed,dataAfterFit,folderName);
save('dataAfterFit','dataAfterFit')

%Calculate all Log2FC's
[dataAfterFit] = calcLog2FC(dataAfterFit);
save('dataAfterFit','dataAfterFit')

%Make heatmap of log fold changes averaged
heatmapFromLog2FCs(dataAfterFit,folderName);

%Make bar plot of raw EC50s with stderr bars, and a dotted horizontal line
%where the parental EC50 is
barPlotsEC50s(dataAfterFit,folderName);

%Make DRC plots for each drug with all cell lines
drugDRCPlots(dataAfterFit,folderName);

%Plot any cell line/drug pair on the same plot (including parentals)


%%
%Everything below here is just tests and stuff

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