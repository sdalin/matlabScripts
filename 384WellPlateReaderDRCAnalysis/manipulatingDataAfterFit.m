% %% Extract all rsquare values into a matrix with the drug_cellLine_date in
% %column 1, and the rsquare in column 2
% 
% rsquaredValues = cell(1,2);
% drugs = fieldnames(dataAfterFit.gofObj.IndFit);
% for drug = 1:size(drugs)
%     cellLines = fieldnames(dataAfterFit.gofObj.IndFit.(drugs{drug}));
%     for cellLine = 1:size(cellLines)
%         rsquaredValues(end+1,1) = {sprintf('%s_%s',drugs{drug},cellLines{cellLine})};
%         rsquaredValues(end,2) = {dataAfterFit.gofObj.IndFit.(drugs{drug}).(cellLines{cellLine}).rsquare};
%     end
% end
% rsquaredValues = rsquaredValues(2:end,1:2);
% 
% %% Make violin plot of r2 values
% clf
% matrsquaredValues = cell2mat(rsquaredValues(:,2));
% for fit = 1:length(matrsquaredValues)
%     nearbyFits = matrsquaredValues(fit)-0.0001 < matrsquaredValues < matrsquaredValues(fit)+0.0001;
%     maxJitter = sum(nearbyFits)/(3*length(matrsquaredValues));
%     jitter = (maxJitter*randn)/3;
%     plot(jitter, matrsquaredValues(fit), 'o', 'MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',3)
%     hold on
% end
% 
% 
% %% Get all EC50s into one big matrix
% drugs = fieldnames(dataAfterFit.rawData);
% drugs = cellfun(@(x) x(6:end),drugs,'UniformOutput',false);
% 
% cellLines = fieldnames(dataAfterFit.rawData.(sprintf('drug_%s',drugs{1})));
% 
% for cellLine = 1:size(cellLines,1)
%         for drug = 1:size(drugs,1)
%             %Find the row/column corresponding to this drug/cell line in
%             %the matrix
%             heatmapRow = find(strcmp(cellLines,cellLines{cellLine}));
%             heatmapColumn = find(strcmp(drugs,(drugs{drug})));
%             
%             %Average of the measured Log2FC's for this cell line/drug pair
%             %First find the column of this cellLine
%             cellLineColumn = find(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC(1,:),cellLines{cellLine}));
%             averageLog2FC = nanmean(cell2mat(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC(2:end,cellLineColumn)));
%             
%             %Put this averaged value into the large matrix in the right
%             %row/column
%             log2FCData(heatmapRow,heatmapColumn) = averageLog2FC;
%         end
% end

%% Get raw data into prism-like output
clear prism
clear tempVector

drugs = fieldnames(dataAfterFit.rawData);
prism.FC = cell(0);

for drug = 1:size(drugs,1)
    cellLines = dataAfterFit.fitParams.(drugs{drug}).Log2FC(1,2:end)';
    prism.(drugs{drug}).cellLineOrder = cellLines;
    %drug doses go in the first column
    prism.(drugs{drug}).doses(:,1) = dataAfterFit.rawData.(drugs{drug}).(cellLines{1})(1:end,1);
    for cellLine = 1:size(cellLines,1)   
        if isfield(dataAfterFit.rawData.(drugs{drug}),cellLines{cellLine})
            %cell line viabilities go in the next columns
            tempToEnterPrism = dataAfterFit.rawData.(drugs{drug}).(cellLines{cellLine})(1:end,2);
            extraLength = size(prism.(drugs{drug}).doses(:,cellLine),1) - size(tempToEnterPrism,1);
            filler = nan(extraLength,1);
            if ~isempty(filler)
                tempToEnterPrism = [tempToEnterPrism;num2cell(filler)];
            end
            
            if size(tempToEnterPrism,1) > size(prism.(drugs{drug}).doses(:,cellLine),1)
                missingLines = size(tempToEnterPrism,1)-size(prism.(drugs{drug}).doses(:,cellLine),1);
                filler = num2cell(nan(missingLines,size(prism.(drugs{drug}).doses,2)-1));
                missingSets = missingLines/16;
                padding = [repmat(dataAfterFit.rawData.(drugs{drug}).(cellLines{1})(1:16,1),missingSets,1),filler];
                prism.(drugs{drug}).doses = [prism.(drugs{drug}).doses;padding];
            end
            prism.(drugs{drug}).doses(:,cellLine+1) = tempToEnterPrism;
            
            
            %if isempty(prism.FC)
            %    prism.FC(1,2:4) = {drugs(drug),NaN,NaN};
            %    prism.FC(2,1) = {cellLines(cellLine)};
            %else
            %    prism.FC(1,end+1:end+3) = {drugs(drug),NaN,NaN};
            %end
            
            cellLineColumn = find(strcmp(cellLines{cellLine},dataAfterFit.fitParams.(drugs{drug}).Log2FC(1,:)));
            tempVector(cellLine+1,2:4) = reshape(dataAfterFit.fitParams.(drugs{drug}).Log2FC(2:4,cellLineColumn),[1,3]);

        end
    end
    tempVector(1,2:4) = {drugs(drug),NaN,NaN};
    
    if isempty(prism.FC)
        prism.FC = tempVector;
    else
        prism.FC(1:6,end+1:end+3) = tempVector(:,2:end);    
    end
    
end
prism.FC(2:6,1) = cellLines;

% 
% clear newPrism
% drugs = fieldnames(prism);
% for drug = 1:size(drugs)
%     newPrism.(drugs{drug})(1:16,1) = prism.(drugs{drug}).doses(1:16,1);
%     for column = 2:size(prism.(drugs{drug}).doses,2)
%         currentViability = prism.(drugs{drug}).doses(:,column);
%         if ~(size(currentViability,1) == 48)
%             sizeDiff = 48 - size(currentViability,1);
%             filler = num2cell(nan(sizeDiff,1));
%             currentViability = [currentViability;filler];
%         end
%         newPrism.(drugs{drug}) = [newPrism.(drugs{drug}),reshape(currentViability,[16,3])];
%     end
% end


    
  