%This function will take fit objects from dataAfterFit struct and make it
%into DRC plots, one for each drug.  Each plot will have all the cell lines
%ever measured for that drug.

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
%DRC plots as descriped above in folder.

function replicateDrugDRCPlots(dataAfterFit,drugsCellLinesThisFolder,folder)
    addpath('/home/sdalin/DRC Analysis for Cluster/kakearney-legendflex-pkg-98b988e/legendflex','/home/sdalin/DRC Analysis for Cluster/kakearney-legendflex-pkg-98b988e/setgetpos_V1.2');

    mkdir(sprintf('%s/matlabOutput/ReplicateDRCPlots',folder))

    allDrugs = drugsCellLinesThisFolder(1,:);
    %drugs = cellfun(@(x) x(6:end),allDrugs,'UniformOutput',false);
    drugs = allDrugs;
    drugs(cellfun(@isempty,drugs)) = [];
    
    for drug = 1:size(drugs,2)
        cellLines = drugsCellLinesThisFolder(2,:);
        cellLines(cellfun(@isempty,cellLines)) = [];
        clf
        
        if ~isfield(dataAfterFit.rawData,(sprintf('drug_%s',drugs{drug})))
            continue
        end
        
        %This struct will hold the handles of the plotted objects
        handles = struct;
        
        %List all fits
        if isfield(dataAfterFit.fitObj.IndFit,(sprintf('drug_%s',drugs{drug})))
            allFits = fieldnames(dataAfterFit.fitObj.IndFit.(sprintf('drug_%s',drugs{drug})));
        end
        
        %Clear legendInfo to not carry over names
        legendInfo = {};
        
        for cellLine = 1:size(cellLines,2)
            
            
            if isempty(cellLines{cellLine})
                noCellLine
                break
            elseif ~isempty(find(~cellfun(@isempty,strfind(dataAfterFit.fitParams.allKiller,sprintf('drug_%s_%s',drugs{drug},cellLines{cellLine}))),1))
                %this line checks to see if this drug/cell line combination
                %was an all-killer If so, continues to the next cell
                %line for that drug.
                continue
            elseif ~isempty(find(~cellfun(@isempty,strfind(dataAfterFit.fitParams.nonKiller,sprintf('drug_%s_%s',drugs{drug},cellLines{cellLine}))),1))
                %this line checks to see if this drug/cell line combination
                %was an non-killer If so, continues to the next cell
                %line for that drug.
                continue
            elseif ~isempty(find(~cellfun(@isempty,strfind(dataAfterFit.fitParams.lessThanThreeDataPoints,sprintf('drug_%s_%s',drugs{drug},cellLines{cellLine}))),1))
                %this line checks to see if this drug/cell line combination
                %had less than 3 points If so, continues to the next cell
                %line for that drug.
            elseif ~isempty(find(~cellfun(@isempty,strfind(dataAfterFit.fitParams.badFit,sprintf('drug_%s_%s',drugs{drug},cellLines{cellLine}))),1))
                %this line checks to see if this drug/cell line combination
                %had rsquare less than 0.8.  If so, continues to the next
                %cell line for that drug.
                continue
            end
            
            replicateFitNames = allFits(strncmp(sprintf('%s_d',cellLines{cellLine}),allFits,length(cellLines{cellLine})+2));
            
            
            %This sets up the colors for the cell lines in the plots
            colorvec = hsv(size(replicateFitNames,1));
            
            %This graphics object array will hold the handles of the lines
            %specifically
            lineHandles = gobjects(size(replicateFitNames,1),1);
            
            if ~isfield(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})),(cellLines{cellLine}))
                continue
            end
            
            concs = dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(cellLines{cellLine})(:,1);
            viability = dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(cellLines{cellLine})(:,2);
            
            %Make concentrations and viabilities into matrixes for ease of
            %processing
            concentrations = cell2mat(concs);
            concsMat = cell2mat(concs);
            viabilities = cell2mat(viability);
            
            for replicate = 1:size(replicateFitNames,1)
                
                %This next line will hide the figures which will probably speeds up
                %running, if you want that
                %p = figure('visible','off')
                fittedHill = dataAfterFit.fitObj.IndFit.(sprintf('drug_%s',drugs{drug})).(replicateFitNames{replicate});
                p = plot(fittedHill,concentrations,viabilities);
                
                %Change colors of curves to see what different replicates are
                %doing
                set(p,'Color',colorvec(replicate,:));
                set(p,'LineWidth',1);
                set(p,'DisplayName',(replicateFitNames{replicate})); 
                
                %Add Rsquare value to the plot
                text(10e-5,0,sprintf('Combined rsquare is %f',dataAfterFit.gofObj.allDataFit.(sprintf('drug_%s',drugs{drug})).(cellLines{cellLine}).rsquare));
                
                %Extract data for making seperate legend file
                legendInfo(replicate)=cellstr(get(p(2),'DisplayName'));
                
                %Save the line handle into a struct and extract the line handle
                %into this nifty graphics objects array
                handles.(replicateFitNames{replicate}) = p;
                lineHandles(replicate,1) = handles.(replicateFitNames{replicate})(2,1);
                
                hold on
            end
  
            axis([-inf,inf,-0.2,1.2])
            set(gca,'XScale','log');
            set(p,'MarkerSize',15)
            xlabel(sprintf('%s Concentration (uM)',drugs{drug}),'FontSize',13);
            ylabel('Viability','FontSize',13);
            legend('off');
            title(sprintf('%s %s',drugs{drug},cellLines{cellLine}),'FontSize',15);
            
            %Make the legends);
            legendCellLines = regexprep(replicateFitNames,'_',' ');
            legendflex(lineHandles(:),legendCellLines,'ncol',1,'anchor',{'ne','ne'},'xscale',0.3,'box','off')
            
            %Save into subfolder of raw data (called 'replicateDRCPlots')
            cd(sprintf('%smatlabOutput/ReplicateDRCPlots',folder))
            set(gcf,'paperpositionmode','manual')
            set(gcf,'paperposition',[-0.25 3.0 9 7])  
            %set(gcf,'Renderer','OpenGL')
            print(sprintf('%s %s replicates DRC plot',drugs{drug},cellLines{cellLine}),'-dpdf');

            cd(sprintf('%s../',folder))
            
            clf
              
        end 
        

    end
            
    clf

end