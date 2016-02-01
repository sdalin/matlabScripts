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
drugNames = {'Doxorubicin','Vincristine','Paclitaxel','Cisplatin'};


dataDir = '/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/Creating Resistant Cells/Round 3/Mondays/CSV_Files/';

numWeeks = length(dir(dataDir))-3;
dose(1,:,1) = [5 6.1 12.2 1700]; %week 1 doses
dose(1,:,2) = [5 6.5 13.5 3500]; %week 2 doses
dose(1,:,3) = [6 7 13.5 3500];%week 3 doses

dateLabel = yyyymmdd(datetime);


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
resWeeks = 3;

for drug = 1:4
    dummy = sum(squeeze(cellLineRankings(:,drug,:))<50,2);
    temporary = find(dummy == resWeeks);
    
    if length(temporary) > length(below50Lines)
        numExtraRows = length(temporary) - length(below50Lines);
        below50Lines = [below50Lines;nan(numExtraRows,drug-1)];
    elseif length(temporary) < length(below50Lines)
        numExtraRows = length(below50Lines) - length(temporary);
        temporary = [temporary;nan(numExtraRows,1)];
    end
    
    below50Lines = [below50Lines, temporary];
end

%Make a matrix with 4 columns, specifying rows of each drug which have
%decreasing rank for the previous resWeeks weeks.

cellLineDiffs = nan(96,1,2);
decreasingLines = [];

for drug = 1:4
    cellLineDiffs = diff(cellLineRankings(:,drug,:),1,3);
    dummy = sum(squeeze(cellLineDiffs)>0);
    temporary = find(dummy == resWeeks-1);
    
    if length(temporary) > length(decreasingLines)
        numExtraRows = length(temporary) - length(decreasingLines);
        below50Lines = [decreasingLines;nan(numExtraRows,drug-1)];
    elseif length(temporary) < length(decreasingLines)
        numExtraRows = length(decreasingLines) - length(temporary);
        temporary = [temporary;nan(numExtraRows,1)];
    end
    
    decreasingLines = [decreasingLines, temporary];
end
    

%%
%Plot rank or z-score vs. week, with selected cell lines from above
%highlighted/colored against grey lines of other cell lines.  Dose can be
%plotted on the right axis as well.

positionVectors= [0.1 0.62 0.3 0.3; 0.6 0.62 0.3 0.3; 0.1 0.1 0.3 0.3; 0.6 0.1 0.3 0.3];

greyco = zeros(7,3);
greyco(:) = 0.8;

normalco = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

%First plot rank vs. week
figure(1)
clf
x = 1:numWeeks;

for drug = 1:4
    subplot(2,2,drug,'Position',positionVectors(drug,:))

    below50LinesNoNAN = below50Lines(~isnan(below50Lines(:,drug)),drug);
    
    y1 = squeeze(cellLineRankings(:,drug,:));
    y2 = squeeze(dose(1,drug,:));
    
    set(gca,'ColorOrder',greyco);
    [hAx,hLine1,hLine2] = plotyy(x,y1,x,y2);
    
    hold on
    
    set(gca,'ColorOrder',normalco)
    y3 = squeeze(cellLineRankings(below50LinesNoNAN,drug,:));    
    sAx = plot(x,y3,'.-','LineWidth',2,'MarkerSize',15);
      
    xlabel('Week')
    ylabel(hAx(1),'Rank') %left y-axis
    formatSpec = '[%s] (nM)';
    ylabel(hAx(2),sprintf(formatSpec,drugNames{drug})) %right y-axis

    axis(hAx(2),[0 numWeeks+1 0 ceil(max(y2))+1])

    hAx(1).YColor = [0 0 0];
    hAx(2).YColor = [1.0000    0.2000    0.2000];
    hAx(1).YLim = [0 100];
    hAx(1).XLim = [0 numWeeks + 1];
    hAx(2).YLim = [0 ceil(max(y2)+(max(y2)/10))];
    hAx(2).XLim = [0 numWeeks + 1];
    hAx(1).YTick = [0 20 40 60 80 100];
    hAx(2).YTick = linspace(0,ceil(max(y2)+(max(y2)/10)),6);

    hLine2.LineStyle = '--';
    hLine2.LineWidth = 2;
    hLine2.Color = [0 0 0];
    hLine2.Marker = 'o';
    hLine2.MarkerFaceColor = [1.0000    0.2000    0.2000];
    hLine2.MarkerEdgeColor = [0 0 0];
    hLine2.MarkerSize = 6;
    
    title(sprintf(drugNames{drug}));

    hold on
end
formatSpec = '%s../%d_Ranks';
print(sprintf(formatSpec,dataDir,dateLabel),'-dpdf');


%Now plot z-score vs. week
figure(2)
clf
x = 1:numWeeks;

for drug = 1:4
    subplot(2,2,drug,'Position',positionVectors(drug,:))

    below50LinesNoNAN = below50Lines(~isnan(below50Lines(:,drug)),drug);
    
    y1 = squeeze(normalizedZscores(:,drug,:));
    y2 = squeeze(dose(1,drug,:));
    
    set(gca,'ColorOrder',greyco);
    [hAx,hLine1,hLine2] = plotyy(x,y1,x,y2);
    
    hold on
    
    set(gca,'ColorOrder',normalco)
    y3 = squeeze(normalizedZscores(below50LinesNoNAN,drug,:));    
    sAx = plot(x,y3,'.-','LineWidth',2,'MarkerSize',15);
      
    xlabel('Week')
    ylabel(hAx(1),'z-score') %left y-axis
    formatSpec = '[%s] (nM)';
    ylabel(hAx(2),sprintf(formatSpec,drugNames{drug})) %right y-axis

    axis(hAx(2),[0 numWeeks+1 0 ceil(max(y2))+1])

    hAx(1).YColor = [0 0 0];
    hAx(2).YColor = [1.0000    0.2000    0.2000];
    
    zscoreLimit = max(abs(min(floor(min(y1)-abs((min(y1)/10))))),max(ceil(max(y1)+(max(y1)/10))));
    
    hAx(1).YLim = [-zscoreLimit zscoreLimit];
    hAx(1).XLim = [0 numWeeks + 1];
    hAx(2).YLim = [0 ceil(max(y2)+(max(y2)/10))];
    hAx(2).XLim = [0 numWeeks + 1];
    hAx(1).YTick = linspace(-zscoreLimit,zscoreLimit,6);
    hAx(2).YTick = linspace(0,ceil(max(y2)+(max(y2)/10)),6);

    hLine2.LineStyle = '--';
    hLine2.LineWidth = 2;
    hLine2.Color = [0 0 0];
    hLine2.Marker = 'o';
    hLine2.MarkerFaceColor = [1.0000    0.2000    0.2000];
    hLine2.MarkerEdgeColor = [0 0 0];
    hLine2.MarkerSize = 6;
    
    title(sprintf(drugNames{drug}));

    hold on
end
formatSpec = '%s../%d_Z-scores';
print(sprintf(formatSpec,dataDir,dateLabel),'-dpdf');


%First plot rank vs. week
figure(3)
clf
x = 1:numWeeks;

for drug = 1:4
    subplot(2,2,drug,'Position',positionVectors(drug,:))

    decreasingLinesNoNAN = decreasingLines(~isnan(decreasingLines(:,drug)),drug);
    
    y1 = squeeze(cellLineRankings(:,drug,:));
    y2 = squeeze(dose(1,drug,:));
    
    set(gca,'ColorOrder',greyco);
    [hAx,hLine1,hLine2] = plotyy(x,y1,x,y2);
    
    hold on
    
    set(gca,'ColorOrder',normalco)
    y3 = squeeze(cellLineRankings(decreasingLinesNoNAN,drug,:));    
    sAx = plot(x,y3,'.-','LineWidth',2,'MarkerSize',15);
      
    xlabel('Week')
    ylabel(hAx(1),'Rank') %left y-axis
    formatSpec = '[%s] (nM)';
    ylabel(hAx(2),sprintf(formatSpec,drugNames{drug})) %right y-axis

    axis(hAx(2),[0 numWeeks+1 0 ceil(max(y2))+1])

    hAx(1).YColor = [0 0 0];
    hAx(2).YColor = [1.0000    0.2000    0.2000];
    hAx(1).YLim = [0 100];
    hAx(1).XLim = [0 numWeeks + 1];
    hAx(2).YLim = [0 ceil(max(y2)+(max(y2)/10))];
    hAx(2).XLim = [0 numWeeks + 1];
    hAx(1).YTick = [0 20 40 60 80 100];
    hAx(2).YTick = linspace(0,ceil(max(y2)+(max(y2)/10)),6);

    hLine2.LineStyle = '--';
    hLine2.LineWidth = 2;
    hLine2.Color = [0 0 0];
    hLine2.Marker = 'o';
    hLine2.MarkerFaceColor = [1.0000    0.2000    0.2000];
    hLine2.MarkerEdgeColor = [0 0 0];
    hLine2.MarkerSize = 6;
    
    title(sprintf(drugNames{drug}));

    hold on
end
formatSpec = '%s../%d_Ranks';
print(sprintf(formatSpec,dataDir,dateLabel),'-dpdf');