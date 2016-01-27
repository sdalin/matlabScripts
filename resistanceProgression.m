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
drugName1 = 'Doxorubicin';
drugName2 = 'Vincristine';
drugName3 = 'Paclitaxel';
drugName4 = 'Cisplatin';

dataDir = '/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/Creating Resistant Cells/Round 3/Mondays/';

numWeeks = length(dir(dataDir))-2;

%%
%First import the PI values for each week.  Rows are different wells,
%columns are different weeks. One matrix per drug.
files = dir(dataDir);

PIData = nan(384,4,numWeeks);

for file = 3:length(files)
    filename = strcat(dataDir,files(file).name);
    PIData(:,:,file-2) = importPIdata(filename);

%     Drug1PI(:,file-2) = PIData(:,1,file-2);
%     Drug2PI(:,file-2) = PIData(:,2,file-2);
%     Drug3PI(:,file-2) = PIData(:,3,file-2);
%     Drug4PI(:,file-2) = PIData(:,4,file-2);
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

% nanRows = any(isnan(Drug1PI),2);
% Drug1PI(nanRows,:) = [];
% nanRows = any(isnan(Drug2PI),2);
% Drug2PI(nanRows,:) = [];
% nanRows = any(isnan(Drug3PI),2);
% Drug3PI(nanRows,:) = [];
% nanRows = any(isnan(Drug4PI),2);
% Drug4PI(nanRows,:) = [];

%%
%Now make similar matrices but with each well normalized like a z-score
%compared to all other wells that week.

normalizedZscores = nan(96,4,numWeeks);

for week = 1:numWeeks
    for drug = 1:4
        normalizedZscores(:,drug,week) = (mean(PIDataFiltered(:,drug,week))-PIDataFiltered(:,drug,week))/std(PIDataFiltered(:,drug,week));
%     Drug1normalized(:,week) = (mean(Drug1PI(:,week))-Drug1PI(:,week))/std(Drug1PI(:,week));
%     Drug2normalized(:,week) = (mean(Drug2PI(:,week))-Drug2PI(:,week))/std(Drug2PI(:,week));
%     Drug3normalized(:,week) = (mean(Drug3PI(:,week))-Drug3PI(:,week))/std(Drug3PI(:,week));
%     Drug4normalized(:,week) = (mean(Drug4PI(:,week))-Drug4PI(:,week))/std(Drug4PI(:,week));
    end
end

%%
%Now make similar matrices but with each well ranked reletive to all other
%wells that week, where the lowest rank is the highest PI value.

cellLineIndices = nan(96,4,numWeeks);
cellLineRankings = nan(96,4,numWeeks);

for week = 1:numWeeks
    for drug = 1:4
        [dummy,cellLineIndices(:,drug,week)] = sort(PIDataFiltered(:,drug,week),'descend');
        cellLineRankings(cellLineIndices(:,drug,week),drug,week) = find(cellLineIndices(:,drug,week));
    end
    
%     [dummy, Drug1indices(:,week)] = sort(Drug1PI(:,week),'descend');
%     Drug1ranked(Drug1indices(:,week),week) = find(Drug1indices(:,week));
%     [dummy, Drug2indices(:,week)] = sort(Drug2PI(:,week),'descend');
%     Drug2ranked(Drug2indices(:,week),week) = find(Drug2indices(:,week));
%     [dummy, Drug3indices(:,week)] = sort(Drug3PI(:,week),'descend');
%     Drug3ranked(Drug3indices(:,week),week) = find(Drug3indices(:,week));
%     [dummy, Drug4indices(:,week)] = sort(Drug4PI(:,week),'descend');
%     Drug4ranked(Drug4indices(:,week),week) = find(Drug4indices(:,week));
end
