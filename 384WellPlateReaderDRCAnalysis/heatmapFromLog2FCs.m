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

function heatmapFromLog2FCs(dataAfterFit,folder)
    
    %First get all the Log2FCs averaged (ignoring []'s/Nan's) and put into
    %a matrix with rows as cell lines and columns as drugs
    drugs = fieldnames(dataAfterFit.rawData);
    drugs = cellfun(@(x) x(6:end),drugs,'UniformOutput',false);

    cellLines = fieldnames(dataAfterFit.rawData.(sprintf('drug_%s',drugs{1})));
   
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
    
    %strip nans because they make clustergram error
    %heatmapMatrix(isnan(heatmapMatrix)) = 0;
    
    %Get values for heatmap scale
    log2ScaleMax = max(max(heatmapMatrix));
    log2ScaleMin = min(min(heatmapMatrix));
    log2Scale = max(log2ScaleMax,abs(log2ScaleMin));
    
    %Make heatmap!
    cellLines = regexprep(cellLines,'_',' ');
    foldChangeHeatmap = clustergram(heatmapMatrix,'RowLabels',cellLines,'ColumnLabels',drugs,'DisplayRange',log2Scale,'Symmetric','true','Colormap',redbluecmap,'RowLabelsRotate',0,'ImputeFun',@knnimpute);
    

    plot(foldChangeHeatmap);
    
    %Somehow this bit is supposed to make the colorbar but its failing
    %miserably. 
    %colormap(redbluecmap(256));
    %imagesc(heatmapMatrix,[-log2Scale log2Scale]);
    %colorbar;
    
    %Save into folder with raw data
    cd(sprintf('%s/matlabOutput',folder))
    set(gcf,'paperpositionmode','auto')
    %set(gcf,'Renderer','OpenGL')
    print('heatmapAllData','-dpdf','-bestfit');
    cd('/Users/sdalin/Dropbox (MIT)/Biology PhD/Matlab Scripts/384WellPlateReaderDRCAnalysis')
    
    close all hidden

end