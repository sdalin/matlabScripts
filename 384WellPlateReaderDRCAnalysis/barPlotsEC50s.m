%This function will take EC50 data from dataAfterFit struct and make it
%into a bar charts, one for each drug.  X-axis will contain each cell line,
%Y-axis will be the EC50 data for each of those cell lines.  Each chart
%will have a dashed horizontal line at the level of the parental EC50

%INPUT: 
%dataAfterFit is a structure with a field for the raw data (subfields
%are drug names, then cell line names, then a matrix with concentrations in the 1st 
%column and normalized viability in the 2nd column. Also fits, with both 
%individual fits with the subfields as drugs,and subfields of that with 
%cell lines & dates of measurement.  All data fit just has cell line names.
%Also a field for paramaters of the fits with similar fields 
%and and a field (fitParams) with subfields of drugnames, as well as lists
%of cell lines that have been fit (dataFit) and cell lines that have had Log2FC
%calculated (dataRun). For drugs, subfields are EC50s and Log2FC which each
%have a cell array with rows as different experiment dates and columns as
%different cell line names.
%Also folder is the path to the folder with the bar charts should be stored.

%OUTPUT:
%Bar charts as descriped above in folder.

function barPlotsEC50s(dataAfterFit,folder)
    mkdir(sprintf('%s/matlabOutput/BarPlots',folder))
    drugs = fieldnames(dataAfterFit.rawData);
    drugs = cellfun(@(x) x(6:end),drugs,'UniformOutput',false);
    
    for drug = 1:size(drugs,1)
        clf
        %Find the measurements of the parental EC50
        columnParental = find(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(1,:),'Parental'));
        parentalEC50Measurements = dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(:,columnParental);
        parentalEC50Average = nanmean(cell2mat(parentalEC50Measurements(2:end)));
        currentEC50s = dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(:,2:end);
        currentEC50s = sortrows(currentEC50s',1)';
        cellLines = currentEC50s(1,:)';
        currentEC50s = currentEC50s(2:end,:);
        currentEC50s(cellfun('isempty',currentEC50s)) = {nan};
        currentEC50s = cell2mat(currentEC50s);
        
        %This next bit makes the dot plot
        dmean = nanmean(currentEC50s,1);                                                                             % Mean
        for cellLine = 1:size(cellLines,1)
            cellLineEC50s = currentEC50s(:,cellLine);
            dci(cellLine) = nanstd(cellLineEC50s)*tinv(0.975,size(cellLineEC50s(~isnan(cellLineEC50s)),1)-1);   % Confidence Intervals
            sem(cellLine) = nanstd(cellLineEC50s)/sqrt(size(cellLineEC50s(~isnan(cellLineEC50s)),1));           %SEM
        end
        xt = [1:size(currentEC50s,2)];                                                                          % X-Ticks
        xtd = repmat(xt, size(currentEC50s,1), 1);                                                              % X-Ticks For Data
        sb = [xt'-ones(size(currentEC50s,2),1)*0.1,  xt'+ones(size(currentEC50s,2),1)*0.1];                     % Short Bar X
        lb = [xt'-ones(size(currentEC50s,2),1)*0.2,  xt'+ones(size(currentEC50s,2),1)*0.2];                     % Long Bar X

        figure(1)
        plot(xt, currentEC50s, 'o', 'MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',3)
        hold on
        for k1 = 1:size(currentEC50s,2)
            plot(lb(k1,:), [1,1]*dmean(k1),'-k',...
                 sb(k1,:), [1,1]*(dmean(k1)-sem(k1)),'-k',...
                 sb(k1,:), [1,1]*(dmean(k1)+sem(k1)),'-k', ...
                 [1,1]*xt(k1),[(dmean(k1)-sem(k1)),(dmean(k1)+sem(k1))],'-k')
        end
        
        plot([0,size(cellLines,1)+1],[1,1]*parentalEC50Average,':k')
        
        hold off
        
        cellLineLabels = regexprep(cellLines,'_',' ');

        set(gca, 'XTick', xt, 'XTickLabel', cellLineLabels,'XTickLabelRotation',90)
        xlabel('Cell Lines')
        ylabel('EC50 (uM)')
        title(drugs{drug})
        if range(range(currentEC50s))
            maxY = max(max(currentEC50s)) + 0.1*range(range(currentEC50s));
            minY = min(min(currentEC50s)) - 0.1*range(range(currentEC50s));
        else
            maxY = max(currentEC50s) + 0.1*range(currentEC50s);
            minY = min(currentEC50s) - 0.1*range(currentEC50s);
        end
        
        axis([0,size(cellLines,1)+1,minY,maxY])
        
        %Save into subfolder of raw data (called 'BarPlots')
        cd(sprintf('%s/matlabOutput/BarPlots',folder))
        set(gcf,'paperpositionmode','auto')
        %set(gcf,'Renderer','OpenGL')
        print(sprintf('%s bar plot',drugs{drug}),'-dpdf','-bestfit');
        cd('/Users/sdalin/Dropbox (MIT)/Biology PhD/Matlab Scripts/384WellPlateReaderDRCAnalysis')

    end
            
    clf

end
