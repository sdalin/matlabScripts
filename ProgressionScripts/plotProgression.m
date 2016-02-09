%%
%Plot rank or z-score vs. week, with selected cell lines from above
%highlighted/colored against grey lines of other cell lines.  Dose can be
%plotted on the right axis as well.

function Plot = plotProgression(numWeeks,dose,drugNames,dataDir,linesToPlot,cellLineRankings,normalizedZscores,selection)

dateLabel = yyyymmdd(datetime);

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
    
    linesToPlotNoNAN = linesToPlot(~isnan(linesToPlot(:,drug)),drug);
    if isempty(linesToPlotNoNAN)
        continue
    end
        
    y1 = squeeze(cellLineRankings(:,drug,:));
    y2 = squeeze(dose(1,drug,:));

    set(groot,'defaultAxesColorOrder',greyco);
    [hAx,hLine1,hLine2] = plotyy(x,y1,x,y2);
    
    hold on
    
    set(gca,'ColorOrder',normalco)
    y3 = squeeze(cellLineRankings(linesToPlotNoNAN,drug,:)); 

    if size(y3,1) == size(y3,2)
        y3 = transpose(y3);
    end
    
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
formatSpec = '%s../%d_%s_Ranks';
print(sprintf(formatSpec,dataDir,dateLabel,selection),'-dpdf');


%Now plot z-score vs. week
figure(2)
clf
x = 1:numWeeks;

for drug = 1:4
    subplot(2,2,drug,'Position',positionVectors(drug,:))

    linesToPlotNoNAN = linesToPlot(~isnan(linesToPlot(:,drug)),drug);
    if isempty(linesToPlotNoNAN)
        continue
    end
    
    y1 = squeeze(normalizedZscores(:,drug,:));
    y2 = squeeze(dose(1,drug,:));
    
    set(groot,'defaultAxesColorOrder',greyco);
    [hAx,hLine1,hLine2] = plotyy(x,y1,x,y2);
    
    hold on
    
    set(gca,'ColorOrder',normalco)
    y3 = squeeze(normalizedZscores(linesToPlotNoNAN,drug,:));    
    if size(y3,1) == size(y3,2)
        y3 = transpose(y3);
    end
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
formatSpec = '%s../%d_%s_Z-scores';
print(sprintf(formatSpec,dataDir,dateLabel,selection),'-dpdf');

Plot = get(gca);
