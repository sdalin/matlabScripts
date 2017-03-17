%Author: Simona Dalin
%Date: 2016/01/19

%This script compiles 384-well PI viability data collected over time into
%plots of which wells have gained resistance to the drug they are treated
%with.  

%Input is the location of a folder containing one csv file for each
%timepoint named as the date in the YYYYMMDD format. In each file there is
%five columns of PI-%'s starting on row 4, One for every well, then one
%each for each of the four drugs used on the 384 well plate FACS'd. (One
%384 well plate is the combination of four 96 well plates). Column 1 is the
%well number, column 2 is the well type which is the drug name in this
%case, column 3 is PI value for all wells, columns 4-7 are PI values for
%each drug with NaN's for wells that don't contain that drug.

%Further in this file you can define 'resWeeks' as the number of weeks you
%want to have plotted cell lines be decreasing in rank or be below rank
%50. You also must define the doses of each drug each week, all these are
%in the first module.

%Outputs are: 4 csv's (one per drug) with the following data: 96-well
%plate ID's, PI value of each well over each week and average killing each week.
%Another similar csv with ranks instead of PI value, where lowest rank is
%highest PI value, this correlates to degree of resistance. Another csv
%with z-scores instead of ranks. Plots of rank and z-score vs. week with
%cell lines with rank below 50 for X weeks highlighted. You can also input
%the dose of each drug at each week and this will be added to the plot.

%%
clear all

%Define some variables
drugNames = {'Doxorubicin','Vincristine','Paclitaxel','Cisplatin'};

%Number of weeks below rank 50 or decreasing in rank that plotted cell
%lines should be.
resWeeks = 10;

dataDir = '/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/Creating Resistant Cells/Round 3/Mondays/CSV_Files/';

dateLabel = yyyymmdd(datetime);

numWeeks = length(dir(dataDir))-3;
dose(1,:,1) = [5 6.1 12.2 1700]; %week 1 doses
dose(1,:,2) = [5 6.5 13.5 3500]; %week 2 doses
dose(1,:,3) = [6 7 13.5 3500];%week 3 doses
dose(1,:,4) = [7.5 10 17.5 3600];%week 4 doses
dose(1,:,5) = [10 12 16 4000];%week 5 doses
dose(1,:,6) = [12.5 10.5 16.5 4150]; %week 6 doses
dose(1,:,7) = [12.5 10.5 16.25 4400]; %week 7 doses
dose(1,:,8) = [15 11 16.5 4750]; %week 8 doses
dose(1,:,9) = [20 11.5 17.5 8500]; %week 9 doses
dose(1,:,10) = [25 12.5 18.6 11000]; %week 10 doses **This week had a 4 day selection rather than 3 days (leitan wedding!)
dose(1,:,11) = [30 15 21 15000]; %week 11 doses (really getting there omg)
dose(1,:,12) = [35 16 21.5 15000]; %week 12 doses (DV from 6/30, PC from 4/15)
dose(1,:,13) = [35 17.5 22.5 15000]; %week 13 doses (DV from 7/8, PC from 5/21)
dose(1,:,14) = [nan nan 24 15250]; %week 14 doses (IQUe was broken for DV week.  PC from 5/27)
dose(1,:,15) = [45 21 26 15250]; %week 15 doses (DV from 7/22, PC from 6/3)
dose(1,:,16) = [nan nan 29 11000]; %week 16 doses (didn't have time to FACS DV because IQue was broken. PC from 6/10)
dose(1,:,17) = [65 25 31.5 15000]; %week 17 doses (DV from 8/5, PC from 6/17)
dose(1,:,18) = [75 27.5 nan nan]; %week 18 doses (DV from 8/21, didn't FACS PC)
%%
%First import the PI values for each week.  Rows are different wells,
%columns are different weeks. One matrix per drug.
files = dir(dataDir);

PIData = nan(384,4,numWeeks);

for file = 4:length(files)
    filename = strcat(dataDir,files(file).name);
    PIData(:,:,file-3) = importPIdata(filename);
end

%%
%Since data is collected in 384 format, there are lots of nan's.  Remove
%those and put everything in 96x4xnumWeeks matrix. 
PIDataFiltered = nan(96,4,numWeeks);

for drug = 1:2
    row96 = 0;
    col96 = 1;
    for group = 0:7
        for row = (drug + (group * 48)): 2 : (drug + 22 + (group * 48));
            PIDataFiltered((col96 + row96),drug,:) = PIData(row,drug,:);
            col96 = col96 + 1;
        end
        col96 = 1;
        row96 = row96 + 12;
    end
end

for drug = 3:4
    row96 = 0;
    col96 = 1;
    for group = 0:7
        for row = (drug + 22 + (group * 48)): 2 : (drug + 44 + (group * 48));
            PIDataFiltered((col96 + row96),drug,:) = PIData(row,drug,:);
            col96 = col96 + 1;
        end
        col96 = 1;
        row96 = row96 +12;
    end
end


%%
%Now make similar matrices but with each well normalized like a z-score
%compared to all other wells that week.

normalizedZscores = nan(96,4,numWeeks);

for week = 1:numWeeks
    for drug = 1:4
        normalizedZscores(:,drug,week) = (PIDataFiltered(:,drug,week)-nanmean(PIDataFiltered(:,drug,week)))/nanstd(PIDataFiltered(:,drug,week));
    end
end

%%
%Now make similar matrices but with each well ranked reletive to all other
%wells that week, where the lowest rank is the highest PI value.

cellLineIndices = PIDataFiltered;
cellLineRankings = PIDataFiltered;

for week = 1:numWeeks
    for drug = 1:4
        loc = find(isnan(PIDataFiltered(:,drug,week)));
        [dummy,cellLineIndices(:,drug,week)] = sort(PIDataFiltered(:,drug,week),'descend');
        cellLineRankings(cellLineIndices(:,drug,week),drug,week) = find(cellLineIndices(:,drug,week));
        cellLineRankings(loc,drug,week) = nan;
        cellLineRankings(:,drug,week) = cellLineRankings(:,drug,week) - length(loc);
    end
end

%%
%Make a matrix with 4 columns, specifying rows of each drug which have rank
%below 50 for the previous resWeeks weeks

%Don't exclude any cell line that has an NaN for a given week withink
%resweeks

below50Lines = [];

for drug = 1:4
    interestingCellLineRankings = cellLineRankings(:,:,numWeeks-resWeeks+1:end);
    %to include cell lines with NaNs, change all Nans to 50 in
    %interestingCellLineRankings
    interestingCellLineRankings50 = interestingCellLineRankings;
    interestingCellLineRankings50(isnan(interestingCellLineRankings)) = 50;
    
    dummy = sum(squeeze(interestingCellLineRankings50(:,drug,:))<=50,2);
    temporary = find(dummy == resWeeks);

    %if ~isempty(below50Lines)
        if size(temporary,1) > size(below50Lines,1)
            numExtraRows = size(temporary,1) - size(below50Lines,1);
            below50Lines = [below50Lines;nan(numExtraRows,drug-1)];
        elseif size(temporary,1) < size(below50Lines,1)
            numExtraRows = size(below50Lines,1) - size(temporary,1);
            temporary = [temporary;nan(numExtraRows,1)];
        end
   % end
    
    below50Lines = [below50Lines, temporary];
end

%Make a matrix with 4 columns, specifying rows of each drug which have
%decreasing rank for the previous resWeeks weeks.

%Don't exclude any cell line that has an NaN for a given week withink
%resweeks

cellLineDiffs = nan(96,1,2);
decreasingLines = [];

for drug = 1:4
    cellLineDiffs = diff(cellLineRankings(:,drug,:),1,3);
    interestingCellLineDiffs = cellLineDiffs(:,:,numWeeks-resWeeks+1:end);
    %to include cell lines with NaNs, change all Nans to 0 in
    %interestingCellLineDiffs
    interestingCellLineDiffs0 = interestingCellLineDiffs;
    interestingCellLineDiffs0(isnan(interestingCellLineDiffs)) = 0;
    
    dummy = sum(squeeze(interestingCellLineDiffs0)<=0,2);    
    temporary = find(dummy == resWeeks-1);
    
    if isempty(temporary) == 1
        temporary = nan;
    end
    
    %if ~isempty(decreasingLines)
        if size(temporary,1) > size(decreasingLines,1)
            numExtraRows = size(temporary,1) - size(decreasingLines,1);
            decreasingLines = [decreasingLines;nan(numExtraRows,drug-1)];
        elseif size(temporary,1) < size(decreasingLines,1)
            numExtraRows = size(decreasingLines,1) - size(temporary,1);
            temporary = [temporary;nan(numExtraRows,1)];           
        end
   % end
    
    decreasingLines = [decreasingLines, temporary];
end
    
%%
%Create lists of which cell lines are below rank 50 or decreasing in rank
% for the previous resWeeks for each drug. Write these into two csv files.

listbelow50Lines = cell(size(below50Lines));
[row,column] = ind2sub([12,8],below50Lines);
asciiCharsCols = char(column+'A'-1);
    for cellLine = 1:size(below50Lines,1)
        for drug = 1:4
            listbelow50Lines{cellLine,drug} = cellstr([asciiCharsCols(cellLine,drug) num2str(row(cellLine,drug))]); 
        end
    end
listbelow50Lines = [drugNames;listbelow50Lines];

    
listDecreasingLines = cell(size(decreasingLines));
[row,column] = ind2sub([12,8],decreasingLines);
asciiCharsCols = char(column+'A'-1);

listDecreasingLines(1,:) = drugNames;
    for cellLine = 1:size(decreasingLines,1)
        for drug = 1:4
            listDecreasingLines{cellLine,drug} = cellstr([asciiCharsCols(cellLine,drug) num2str(row(cellLine,drug))]); 
        end
    end
listDecreasingLines = [drugNames;listDecreasingLines];

%%
%Plot everything first for below50Lines
Plot = plotProgression(numWeeks,dose,drugNames,dataDir,below50Lines,cellLineRankings,normalizedZscores,PIDataFiltered,'below50');

%Plot everything now for decreasingLines
Plot = plotProgression(numWeeks,dose,drugNames,dataDir,decreasingLines,cellLineRankings,normalizedZscores,PIDataFiltered,'decreasing');

