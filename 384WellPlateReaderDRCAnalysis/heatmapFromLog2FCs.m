%This function will take Log2FC data from dataAfterFit struct and make it
%into a heatmap with cell lines on the y-axis, and drugs on the x-axis.  It
%will do it as a clustergram, with heierarchical clustering, to see which
%drugs/cell line are most similar to each other.

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
%Also folder is the path to the folder with the heatmap should be stored.

%OUTPUT:
%A heatmap as descriped above in folder.

function [heatmapMatrixCleaned,cellLines,drugsCleaned] = heatmapFromLog2FCs(dataAfterFit,folder,stars,varargin)
    
    close all hidden
    %First get all the Log2FCs averaged (ignoring []'s/Nan's) and put into
    %a matrix with rows as cell lines and columns as drugs
    drugs = fieldnames(dataAfterFit.rawData);
    drugs = cellfun(@(x) x(6:end),drugs,'UniformOutput',false);

    cellLines = fieldnames(dataAfterFit.rawData.(sprintf('drug_%s',drugs{1})));
    
    %Extract just cell lines resistant to a particular drug, to make
    %individual heatmaps
    if nargin > 3
        cellLinesOneDrug = cellLines(strncmpi(cellLines,varargin{1},3));
        cellLinesDMSO = cellLines(strncmpi(cellLines,'DMSO',4));
        cellLines = vertcat(cellLinesOneDrug,cellLinesDMSO);
    end
   
   
    %Initialize the matrix
    heatmapMatrix = [];
    
    for cellLine = 1:size(cellLines,1)
        for drug = 1:size(drugs,1)
            %Find the row/column corresponding to this drug/cell line in
            %the matrix
            heatmapRow = find(strcmp(cellLines,cellLines{cellLine}));
            heatmapColumn = find(strcmp(drugs,(drugs{drug})));
            
            %Average of the measured Log2FC's for this cell line/drug pair
            %First find the column of this cellLine
            cellLineColumn = find(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC(1,:),cellLines{cellLine}));
            averageLog2FC = nanmean(cell2mat(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC(2:end,cellLineColumn)));
            
            %Put this averaged value into the large matrix in the right
            %row/column
            heatmapMatrix(heatmapRow,heatmapColumn) = averageLog2FC;
        end
    end
    
    %Remove any rows/columns that are all nan, and corresponding row/column
    %label
    emptyRowsColumns = 1;
    while emptyRowsColumns
        drugs(~any(~isnan(heatmapMatrix))') = [];
        heatmapMatrix(:,~any(~isnan(heatmapMatrix))) = [];
        cellLines(~any(~isnan(heatmapMatrix),2)) = [];
        heatmapMatrix(~any(~isnan(heatmapMatrix),2),:) = [];
        
        if ~isempty(find(~any(~isnan(heatmapMatrix)),1)) || ~isempty(find(~any(~isnan(heatmapMatrix),2),1))
            emptyRowsColumns = 1;
        elseif isempty(find(~any(~isnan(heatmapMatrix)),1)) && isempty(find(~any(~isnan(heatmapMatrix),2),1))
            emptyRowsColumns = 0;
        end
    end
    
    
    
    
    %Get values for heatmap scale
    log2ScaleMax = max(max(heatmapMatrix));
    log2ScaleMin = min(min(heatmapMatrix));
    log2Scale = max(log2ScaleMax,abs(log2ScaleMin));
    
    %Make heatmap!
    cellLines = regexprep(cellLines,'_',' ');

    %Remove drugs with data for less than half the cell lines
    if size(heatmapMatrix,1) == 1
        heatmapMatrixCleaned = heatmapMatrix;
        drugsCleaned = drugs;
    else
        heatmapMatrixCleaned = heatmapMatrix(:,sum(~isnan(heatmapMatrix))>size(heatmapMatrix,1)/2);
        drugsCleaned = drugs(sum(~isnan(heatmapMatrix))>size(heatmapMatrix,1)/2);
    end
    
    %Remove drugs with log2(FC) >|1| in >2 DMSO control lines
    heatmapMatrixDMSOLines = heatmapMatrixCleaned(strncmp('DMSO',cellLines,4),:);
    heatmapMatrixNoArtifact = heatmapMatrixCleaned(:,~(sum(heatmapMatrixDMSOLines < -1) > 2));
    drugsNoArtifact = drugsCleaned(~(sum(heatmapMatrixDMSOLines < -1) > 2));
    
    if size(heatmapMatrix,1) == 1
        clustVal = 2;
    else
        clustVal = 3;
    end
        
    foldChangeHeatmap = clustergram(heatmapMatrixNoArtifact,'RowLabels',cellLines,'ColumnLabels',drugsNoArtifact,'DisplayRange',log2Scale,'Symmetric','true','Colormap',redbluecmap,'RowLabelsRotate',0,'ImputeFun',@knnimpute,'Cluster',clustVal);
    
    
    p = plot(foldChangeHeatmap);
    
    %adjust the colorbar to look nice
    colorBar = findobj('Tag','HeatMapColorbar');
    colorBar.FontSize = 10;
    colorBar.Position = [0.075,0.125,0.01,0.6];
    
    hold on
    for cellLine = 1:size(foldChangeHeatmap.RowLabels,1)
        for drug = 1:size(foldChangeHeatmap.ColumnLabels,2)
            [drugRow,col] = find(strcmp(stars,sprintf('%s',foldChangeHeatmap.ColumnLabels{drug})));
            [heatmapDrugRow,col] = find(strcmp(drugsNoArtifact,sprintf('%s',foldChangeHeatmap.ColumnLabels{drug})));
            [row,cellLineColumn] = find(strcmp(stars,sprintf('%s',foldChangeHeatmap.RowLabels{cellLine})));
            [heatmapCellLineColumn,tmp] = find(strcmp(cellLines,sprintf('%s',foldChangeHeatmap.RowLabels{cellLine})));
            
            if ~isempty(stars{drugRow,cellLineColumn}) && ~isnan(heatmapMatrixNoArtifact(heatmapCellLineColumn,heatmapDrugRow))
                if nargin > 3
                    starText = text(p,drug,cellLine,char(stars{drugRow,cellLineColumn}),'Color',[0.5 0.5 0.5],'HorizontalAlignment','center','VerticalAlignment','middle');
                else
                    starText = text(p,drug,cellLine-0.5,char(stars{drugRow,cellLineColumn}),'Color',[0.5 0.5 0.5],'HorizontalAlignment','center','VerticalAlignment','middle');
                end
            end
        end
    end
    
    y = size(foldChangeHeatmap.RowLabels,1)+(size(foldChangeHeatmap.RowLabels,1)*((130-112)/112));
    x = size(foldChangeHeatmap.ColumnLabels,2)+(size(foldChangeHeatmap.ColumnLabels,2)*((-3-20)/20));
    log2FCText = text(p,x,y,'Log_{2}\DeltaEC50','HorizontalAlignment','center','Interpreter','tex','FontSize',15);
    
    y = size(foldChangeHeatmap.RowLabels,1)+(size(foldChangeHeatmap.RowLabels,1)*((115-112)/112));
    x = size(foldChangeHeatmap.ColumnLabels,2)+(size(foldChangeHeatmap.ColumnLabels,2)*((-5.2-20)/20)); 
    colResText = text(p,x,y,'ColRes','HorizontalAlignment','center','Color',[0.4039 0 0.1216]);
    
    y = size(foldChangeHeatmap.RowLabels,1)+(size(foldChangeHeatmap.RowLabels,1)*((-1-112)/112));
    x = size(foldChangeHeatmap.ColumnLabels,2)+(size(foldChangeHeatmap.ColumnLabels,2)*((-5.2-20)/20)); 
    colSenText = text(p,x,y,'ColSen','HorizontalAlignment','center','Color',[0.0196 0.1882 0.3804]);
    
    
    %Save into folder with raw data
    cd(sprintf('%smatlabOutput',folder))
    set(gcf,'paperpositionmode','manual')
    set(gcf,'PaperUnits','inches')
    set(gcf,'paperposition',[-0.30 -0.25 9 12])
    set(gcf,'Renderer','OpenGL')
    if nargin > 3
        print(sprintf('heatmapAllData_%s',varargin{1}),'-dpdf');
    else
        print('heatmapAllData','-dpdf');
    end
    cd(sprintf('%s../',folder))
    
    clf
    close all hidden
    
end
