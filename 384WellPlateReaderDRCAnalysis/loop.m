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
folderName = sprintf('%s/Raw Data/',directory);
list=dir(sprintf('%s*.csv',folderName));

%% The following is for cases where each plate is one csv file and there is a directory linking assay & compound plate barcodes 

%Read in drug names, KI numbers, and targets into bigstructNormed.Info
infoFileName = dir(sprintf('%s/KI-HTS.Selleck.384wellarray.20141007.xlsx',directory));
[num,drugData,raw] = xlsread(infoFileName.name,'Chemical Data','A2:F393');
bigstructNormed.Info = [drugData(:,1:2),drugData(:,5)];

%Read in the cell line/assay plate barcode/compound plate barcode data
directoryFileName = dir(sprintf('%s/Simona_Selleck_20171218.xlsx',directory));
[num,cellLinePlateCompoundBarcodes,raw] = xlsread(directoryFileName.name,'Sheet1','A2:D61');

%Read in each assay plate layout & drug concentrations
%Find number of assay plates
[status,sheets] = xlsfinfo(infoFileName.name);
numOfAssayPlates = numel(sheets) - 1;

for assayPlate = 2:numOfAssayPlates+1
    %Grab the doses
    [num,text,doses] = xlsread(infoFileName.name,sheets{numOfAssayPlates},'B21:Y37');
    
    %Store doses in bigstructNormed.Info
    bigstructNormed.Directory.(sheets{numOfAssayPlates}).doses = doses;
    
    %Grab the drugNames
    [num,drugNames,raw] = xlsread(infoFileName.name,sheets{numOfAssayPlates},'B4:Y19');
    
    %Store drug names in BigstructNormed.Info
    bigstructNormed.Directory.(sheets{numOfAssayPlates}).drugNames = drugNames;
end
    
%%

for i = 1:length(list);
    filename = sprintf('%s%s',folderName,list(i).name);
    bigstructNormed.(list(i).name(1:end-5)) = ReadNormPlateReaderHTSDRCData(filename,cellLinePlateCompoundBarcodes,bigstructNormed.Directory);
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
[drugsCellLinesThisFolder] = hillFitAllHTS(bigstructNormed,folderName);

%Calculate all Log2FC's
load('dataAfterFit','dataAfterFit');
[dataAfterFit] = calcLog2FCHTS(dataAfterFit);
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
load(sprintf('%s/Raw Data/matlabOutput/drugsCellLinesThisFolder',directory));
replicateDrugDRCPlotsHTS(dataAfterFit,drugsCellLinesThisFolder,folderName);



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