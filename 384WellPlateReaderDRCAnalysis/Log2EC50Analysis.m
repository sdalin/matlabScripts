%Find optimal k for kmeans clustering, do PCA, then re-do optimal k, for
%all data and then for each drug-resistance group individually

function [eva,differenceReconData,differenceCoeffsPCAPPCA,differenceCoeffsALSPPCA,differenceCoeffsALSPCA,evaPCA,evaPPCA,evaALS] = Log2EC50Analysis(dataAfterFit,folder,dataType)
    %% Extracting desired data
    drugs = fieldnames(dataAfterFit.rawData);
    drugs = cellfun(@(x) x(6:end),drugs,'UniformOutput',false);
    
    cellLines = fieldnames(dataAfterFit.rawData.(sprintf('drug_%s',drugs{1})));
    
    if ~(strcmp(dataType,'AllData'))
        cellLinesOneDrug = cellLines(strncmpi(cellLines,dataType,3));
        cellLines = cellLinesOneDrug;
    end

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
            
            %Remove drugs with data for less than half the cell lines
            log2FCDataCleaned = log2FCData(:,sum(~isnan(log2FCData))>size(log2FCData,1)/2);
            drugsCleaned = drugs(sum(~isnan(log2FCData))>size(log2FCData,1)/2);
        end
    end
    
    %% Determine optimal number of clusters & create elbow plot
    
    [eva] = kmeansAnalysis(log2FCDataCleaned,folder,dataType);
    
    if eva.OptimalK > 13
        k = eva.OptimalK;
    else
        k = 4;
    end
    
    idxRaw = kmeans(log2FCDataCleaned,k);

    
    %% PCA
    [wcoeff,score,latent,tsquared,explained] = pca(log2FCDataCleaned);
    [wcoeffALS,scoreALS,latentALS,tsquaredALS,explainedALS] = pca(log2FCDataCleaned,'algorithm','als');
    
    
    %% Probabalistic PCA
    
    [wcoeffPPCA,scorePPCA,latentPPCA,mu,v,S] = ppca(log2FCDataCleaned,min(size(log2FCDataCleaned)));
    
    explained = 100*latentPPCA/sum(latentPPCA);
        

    differenceReconData = S.Recon - log2FCDataCleaned;
    differenceCoeffsPCAPPCA = subspace(wcoeff,wcoeffPPCA);
    differenceCoeffsALSPPCA = subspace(wcoeffALS,wcoeffPPCA);
    differenceCoeffsALSPCA = subspace(wcoeffALS,wcoeff);
    %% Determine optimal number of clusters on dimensionality reduced data
    [evaPCA] = kmeansAnalysis(score(:,[1:3]),folder,dataType,'PCA');
    if evaPCA.OptimalK > 13
        k = evaPCA.OptimalK;
    else
        k = 4;
    end
    idxPCA = kmeans(log2FCDataCleaned,k);
    
    [evaPPCA] = kmeansAnalysis(scorePPCA(:,[1:3]),folder,dataType,'PPCA');
    
    if 1 < evaPPCA.OptimalK < 14
        k = evaPPCA.OptimalK;
    else
        k = 4;
    end
    dataType
    idxPPCA = kmeans(log2FCDataCleaned,k);
    
    
    [evaALS] = kmeansAnalysis(scoreALS(:,[1:3]),folder,dataType,'PCAALS');
    
    if evaALS.OptimalK > 13
        k = evaALS.OptimalK;
    else
        k = 4;
    end
    
    idxALS = kmeans(log2FCDataCleaned,k);
    
    
    %% Make biplots with points labeled according to kmeans on raw data, kmeans on PCA, kmeans on PPCA, kemans on PCAALS
    makeParetoandBiPlots(explained,dataType,'PCA-Raw',folder,wcoeffPPCA,scorePPCA,drugsCleaned,cellLines,idxRaw)
    makeParetoandBiPlots(explained,dataType,'PCA',folder,wcoeffPPCA,scorePPCA,drugsCleaned,cellLines,idxPCA)
    makeParetoandBiPlots(explained,dataType,'PPCA',folder,wcoeffPPCA,scorePPCA,drugsCleaned,cellLines,idxPPCA)
    makeParetoandBiPlots(explained,dataType,'ALS',folder,wcoeffPPCA,scorePPCA,drugsCleaned,cellLines,idxALS)


end