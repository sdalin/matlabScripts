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
directory = pwd;
bigstructNormed = struct;
%this is going to be a nested struct containing ALL the data.  First level
%is different days of data collection.  Second level is separate plates on each day and one plate with 'info' ie, drug names and concentrations.
folderName = sprintf('%s/Raw Data weeks of 161118 170213 170227 170313/',directory);
list=dir(sprintf('%s*.xlsx',folderName));
for i = 1:length(list);
    filename = sprintf('%s%s',folderName,list(i).name);
    bigstructNormed.(sprintf('date%s',filename(end-14:end-9))) = ReadNormPlateReaderDRCData(filename);
end

%Get fits of all raw data in this folder  
%First make some folders to store images in
if ~isdir(sprintf('%s/matlabOutput',folderName))
    mkdir(sprintf('%s/matlabOutput',folderName));
end

if ~isdir(sprintf('%s/matlabOutput/rawDRCs',folderName))
    mkdir(sprintf('%s/matlabOutput/rawDRCs',folderName));
end

if ~isdir(sprintf('%s/matlabOutput/rejectedFits',folderName))
    mkdir(sprintf('%s/matlabOutput/rejectedFits',folderName));
end

%folderName = sprintf('%s..',folderName);

%Calculate fits of individual data, all replicates, and put graphs of individual fits in the
%folder above.
[drugsCellLinesThisFolder] = hillFitAll(bigstructNormed,folderName);

%Calculate all Log2FC's
load('dataAfterFit','dataAfterFit');
[dataAfterFit] = calcLog2FC(dataAfterFit);
save('dataAfterFit','dataAfterFit')

%Make bar plot of raw EC50s with stderr bars, and a dotted horizontal line
%where the parental EC50 is
[h,p,stars] = barPlotsEC50s(dataAfterFit,folderName);

%Make heatmap of log fold changes averaged
[heatmapMatrixAll,cellLinesAll,drugsCleanedAll] = heatmapFromLog2FCs(dataAfterFit,folderName,stars);
[heatmapMatrixDox,cellLinesDox,drugsCleanedDox] = heatmapFromLog2FCs(dataAfterFit,folderName,stars,'Dox');
[heatmapMatrixVin,cellLinesVin,drugsCleanedVin] = heatmapFromLog2FCs(dataAfterFit,folderName,stars,'Vin');
[heatmapMatrixPac,cellLinesPac,drugsCleanedPac] = heatmapFromLog2FCs(dataAfterFit,folderName,stars,'Pac');
[heatmapMatrixCis,cellLinesCis,drugsCleanedCis] = heatmapFromLog2FCs(dataAfterFit,folderName,stars,'Cis');

%Make DRC plots for each drug with all cell lines in this particular folder
load(sprintf('%s/Raw Data weeks of 161118 170213 170227 170313/matlabOutput/drugsCellLinesThisFolder',directory));
replicateDrugDRCPlots(dataAfterFit,drugsCellLinesThisFolder,folderName);



%% This bit can't be run on the cluster
%Select which cell lines and drugs to plot
[selectedDrugs,selectedCellLines] = selectDataToPlot(dataAfterFit);

%Make plots with all replicates of selected cell line/drug pair on the same
%plot, and save it into a folder called replicateDRCPlots
drugDRCPlots(dataAfterFit,selectedDrugs,selectedCellLines,folderName);



%% I'm not sure what this bit is.
% %Everything below here is just tests and stuff
% 
% %Get values for heatmap scale
% log2ScaleMax = max(max(log2FoldChangesnoEu));
% log2ScaleMin = min(min(log2FoldChangesnoEu));
% log2Scale = max(log2ScaleMax,abs(log2ScaleMin));
% 
% %Make heatmap!
% foldChangeHeatmap = clustergram(log2FoldChangesnoEu','RowLabels',drugNames,'ColumnLabels',cellLinesnoEu,'DisplayRange',log2Scale,'Symmetric','true','Colormap',redbluecmap,'ColumnLabelsRotate',0);
% plot(foldChangeHeatmap);
% colormap(redbluecmap(256));
% imagesc(log2FoldChangesnoEu',[-log2Scale log2Scale]);
% colorbar;