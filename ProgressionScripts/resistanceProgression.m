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
%want to have plotted cell lines be decreasing in rank or be below rank 50.

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
resWeeks = 3;

dataDir = '/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/Creating Resistant Cells/Round 3/Mondays/CSV_Files/';

numWeeks = length(dir(dataDir))-3;
dose(1,:,1) = [5 6.1 12.2 1700]; %week 1 doses
dose(1,:,2) = [5 6.5 13.5 3500]; %week 2 doses
dose(1,:,3) = [6 7 13.5 3500];%week 3 doses


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

below50Lines = [];

for drug = 1:4
    interestingCellLineRankings = cellLineRankings(:,:,numWeeks-resWeeks+1:end);
    dummy = sum(squeeze(interestingCellLineRankings(:,drug,:))<50,2);
    temporary = find(dummy == resWeeks);

    if ~isempty(below50Lines)
        if length(temporary) > length(below50Lines)
            numExtraRows = length(temporary) - length(below50Lines);
            below50Lines = [below50Lines;nan(numExtraRows,drug-1)];
        elseif length(temporary) < length(below50Lines)
            numExtraRows = length(below50Lines) - length(temporary);
            temporary = [temporary;nan(numExtraRows,1)];
        end
    end
    
    below50Lines = [below50Lines, temporary];
end

%Make a matrix with 4 columns, specifying rows of each drug which have
%decreasing rank for the previous resWeeks weeks.

cellLineDiffs = nan(96,1,2);
decreasingLines = [];

for drug = 1:4
    cellLineDiffs = diff(cellLineRankings(:,drug,:),1,3);
    interestingCellLineDiffs = cellLineDiffs(:,:,numWeeks-resWeeks+1:end);
    dummy = sum(squeeze(interestingCellLineDiffs)<0,2);    
    temporary = find(dummy == resWeeks-1);
    
    if ~isempty(decreasingLines)
        if length(temporary) > length(decreasingLines)
            numExtraRows = length(temporary) - length(decreasingLines);
            decreasingLines = [decreasingLines;nan(numExtraRows,drug-1)];
        elseif length(temporary) < length(decreasingLines)
            numExtraRows = length(decreasingLines) - length(temporary);
            temporary = [temporary;nan(numExtraRows,1)];
        end
    end
    
    decreasingLines = [decreasingLines, temporary];
end
    

%%
%Plot everything first for below50Lines
Plot = plotProgression(numWeeks,dose,drugNames,dataDir,below50Lines,cellLineRankings,normalizedZscores,'below50');

%Plot everything now for decreasingLines
Plot = plotProgression(numWeeks,dose,drugNames,dataDir,decreasingLines,cellLineRankings,normalizedZscores,'decreasing');
