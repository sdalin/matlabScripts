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

function [h,p,stars] = barPlotsEC50s(dataAfterFit,folder)
    allCellLines = fieldnames(dataAfterFit.rawData.('drug_Doxorubicin'));
    
    figure('visible','off');
    
    mkdir(sprintf('%s/matlabOutput/BarPlots',folder))
    drugs = fieldnames(dataAfterFit.rawData);
    drugs = cellfun(@(x) x(6:end),drugs,'UniformOutput',false);
    
    %initialize cell array to store stats
    h = cell(size(drugs,1)+1,size(allCellLines,1)+1);
    h(:,1) = ['NaN';drugs];
    h(1,:) = ['NaN';allCellLines]';
    
    p = h;
    stars = h;
    
    for drug = 1:size(drugs,1)
        close all hidden
        hold off
        %Find the measurements of the parental EC50
        columnParental = find(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(1,:),'Parental'));
        parentalEC50Measurements = dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(:,columnParental);
        parentalEC50Average = nanmean(cell2mat(parentalEC50Measurements(2:end)));
        currentEC50s = dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(:,2:end);
        %currentEC50s = sortrows(currentEC50s',1)';
        cellLines = currentEC50s(1,:)';
        currentEC50s = currentEC50s(2:end,:);
        currentEC50s(cellfun('isempty',currentEC50s)) = {nan};
        currentEC50s = cell2mat(currentEC50s);
        
        %2 sample t-test with bonferroni correction of each cell line's
        %EC50s vs. parental EC50s
        for cellLine = 1:size(currentEC50s,2)
            if ~sum(~isnan(currentEC50s(:,cellLine)))
                continue
            end
            
            newAlpha05 = 0.05/(sum(~isnan(currentEC50s(:,cellLine))));
            newAlpha01 = 0.01/(sum(~isnan(currentEC50s(:,cellLine))));
            newAlpha001 = 0.001/(sum(~isnan(currentEC50s(:,cellLine))));
            
            newAlpha = [newAlpha05,newAlpha01,newAlpha001];
            
            
            [drugRow,col] = find(strcmp(h,sprintf('%s',drugs{drug})));
            [row,cellLineColumn] = find(strcmp(h,sprintf('%s',cellLines{cellLine})));
            
            for alpha = 1:size(newAlpha,2)
                [h{drugRow,cellLineColumn},p{drugRow,cellLineColumn}] = ttest2(currentEC50s(:,cellLine),cell2mat(parentalEC50Measurements(2:end,:)),'Alpha',newAlpha(alpha));
                
                if strcmp(cellLines{cellLine},'Parental')
                    continue;
                elseif cell2mat(h(drugRow,cellLineColumn))
                    stars{drugRow,cellLineColumn} = (repmat('*',[1,alpha]));
                end
            end
        end
        
        %remove _ from cell line names in stars
        stars(1,:) = cellfun(@(x) strrep(x,'_',' '),stars(1,:),'Uniform',false);
        
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

        set(gcf,'visible','off')
        plot(xt, currentEC50s, 'o', 'MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',3)
        hold on
        for k1 = 1:size(currentEC50s,2)
            plot(lb(k1,:), [1,1]*dmean(k1),'-k',...
                 sb(k1,:), [1,1]*(dmean(k1)-sem(k1)),'-k',...
                 sb(k1,:), [1,1]*(dmean(k1)+sem(k1)),'-k', ...
                 [1,1]*xt(k1),[(dmean(k1)-sem(k1)),(dmean(k1)+sem(k1))],'-k')
        end
        
        plot([0,size(cellLines,1)+1],[1,1]*parentalEC50Average,':k')
        
        
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
        set(gcf,'paperpositionmode','manual')
        set(gcf,'paperposition',[-0.75 0.1 10 11])
        %set(gcf,'Renderer','OpenGL')
        print(sprintf('%s bar plot',drugs{drug}),'-dpdf');
        cd(sprintf('%s../',folder))

        hold off

    end
            
    close all hidden
    
    
    
end
