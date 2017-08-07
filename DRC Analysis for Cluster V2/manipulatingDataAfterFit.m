%Extract all rsquare values into a matrix with the drug_cellLine_date in
%column 1, and the rsquare in column 2

rsquaredValues = cell(1,2);
drugs = fieldnames(dataAfterFit.gofObj.IndFit);
for drug = 1:size(drugs)
    cellLines = fieldnames(dataAfterFit.gofObj.IndFit.(drugs{drug}));
    for cellLine = 1:size(cellLines)
        rsquaredValues(end+1,1) = {sprintf('%s_%s',drugs{drug},cellLines{cellLine})};
        rsquaredValues(end,2) = {dataAfterFit.gofObj.IndFit.(drugs{drug}).(cellLines{cellLine}).rsquare};
    end
end
rsquaredValues = rsquaredValues(2:end,1:2);

%% Make violin plot of r2 values
clf
matrsquaredValues = cell2mat(rsquaredValues(:,2));
for fit = 1:length(matrsquaredValues)
    nearbyFits = matrsquaredValues(fit)-0.0001 < matrsquaredValues < matrsquaredValues(fit)+0.0001;
    maxJitter = sum(nearbyFits)/(3*length(matrsquaredValues));
    jitter = (maxJitter*randn)/3;
    plot(jitter, matrsquaredValues(fit), 'o', 'MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',3)
    hold on
end


%% Get all EC50s into one big matrix
drugs = fieldnames(dataAfterFit.rawData);
drugs = cellfun(@(x) x(6:end),drugs,'UniformOutput',false);

cellLines = fieldnames(dataAfterFit.rawData.(sprintf('drug_%s',drugs{1})));

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
            log2FCData(heatmapRow,heatmapColumn) = averageLog2FC;
        end
    end
