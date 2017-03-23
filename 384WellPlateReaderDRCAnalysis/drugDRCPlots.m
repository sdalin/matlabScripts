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

function drugDRCPlots(dataAfterFit,selectedDrugs,selectedCellLines,folder)
    addpath('/Users/sdalin/Dropbox (MIT)/Biology PhD/Matlab Scripts/kakearney-legendflex-pkg-98b988e/legendflex','/Users/sdalin/Dropbox (MIT)/Biology PhD/Matlab Scripts/kakearney-legendflex-pkg-98b988e/setgetpos_V1.2');

    mkdir(sprintf('%s/matlabOutput/DRCPlots',folder))
    drugs = selectedDrugs;
    drugs = cellfun(@(x) x(6:end),drugs,'UniformOutput',false);
    allDrugs = fieldnames(dataAfterFit.rawData);
    
    for drug = 1:size(drugs,1)
        allCellLines = fieldnames(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})));
        cellLines = intersect(allCellLines,selectedCellLines);
        
        if isempty(cellLines) %If this drug/cell line pair doesn't exist, skip it.
            continue
        end
        
        clf
        %This sets up the colors for the cell lines in the plots
        colorvec = hsv(size(cellLines,1));
        
        %This struct will hold the handles of the plotted objects
        handles = struct;
        
        %This graphics object array will hold the handles of the lines
        %specifically
        lineHandles = gobjects(size(cellLines,1),1);
        
        for cellLine = 1:size(cellLines,1)
            fittedHill = dataAfterFit.fitObj.allDataFit.(sprintf('drug_%s',drugs{drug})).(cellLines{cellLine});
            concs = dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(cellLines{cellLine})(:,1);
            viability = dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(cellLines{cellLine})(:,2);
            
            %Make concentrations and viabilities into matrixes for ease of
            %processing
            concentrations = unique(cell2mat(concs));
            concsMat = cell2mat(concs);
            viabilities = cell2mat(viability);
            
            %Calculate mean and SEM of all viabilities for each
            %concentration
            sem = nan(size(concentrations,1),1);
            means = nan(size(concentrations,1),1);
            for concentration = 1:size(concentrations)
                currentViabilities = viabilities(concsMat == concentrations(concentration));
                sem(concentration,1) = nanstd(currentViabilities)/sqrt(size(currentViabilities(~isnan(currentViabilities)),1));           
                means(concentration,1) = nanmean(currentViabilities);
            end
            
            %This next line will hide the figures which will probably speeds up
            %running, if you want that
            %p = figure('visible','off');
            p = plot(fittedHill,concentrations,means);
       
            %Change the way it looks so its easily readable
            %Set the parental curve to black and other curves to rainbow
            %colors
            if strcmp(cellLines{cellLine},'Parental')
                set(p,'Color',[0,0,0]);
                set(p,'LineWidth',3);
            else
                set(p,'Color',colorvec(cellLine,:));
                set(p,'LineWidth',1);
            end
            
            cellLineLabels = regexprep(cellLines,'_',' ');
            
            set(p,'DisplayName',(cellLineLabels{cellLine}));  
            axis([-inf,inf,-0.2,1.2])
            set(gca,'XScale','log');
            set(p,'MarkerSize',15)
            xlabel(sprintf('%s Concentration (uM)',drugs{drug}),'FontSize',13);
            ylabel('Viability','FontSize',13);
            legend('off');
            title(sprintf(drugs{drug}),'FontSize',15);
            
            %Extract data for making seperate legend file
            legendInfo(cellLine)=cellstr(get(p(2),'DisplayName'));
            
            %Save the line handle into a struct and extract the line handle
            %into this nifty graphics objects array
            handles.(cellLines{cellLine}) = p;
            lineHandles(cellLine,1) = handles.(cellLines{cellLine})(2,1);
            
            %Keep plot here for the error bars
            hold on
            
            %Put in the error bars
             if strcmp(cellLines{cellLine},'Parental')
                errorbar(concentrations,means,sem,'Color',[0,0,0],'Linestyle','none','LineWidth',2)
             else
                errorbar(concentrations,means,sem,'Color',colorvec(cellLine,:),'Linestyle','none')
            end
            
            %keep keeping the plot on for the next plot
            hold on
            
        end 
        
        %Make the legends);
        legendCellLines = regexprep(cellLines,'_',' ');
        legendflex(lineHandles(:),legendCellLines,'ncol',2,'anchor',{'ne','se'},'buffer',[105,-200],'xscale',0.3,'box','off')
        %Nudge everything left a smidge so it fits in the plot
        pos = get(gca, 'Position');
        xoffset = -0.07;
        pos(1) = pos(1) + xoffset;
        set(gca, 'Position', pos)
          
        %Save into subfolder of raw data (called 'DRCPlots')
        cd(sprintf('%s/matlabOutput/DRCPlots',folder))
        set(gcf,'paperpositionmode','auto')
        %set(gcf,'Renderer','OpenGL')
        if sum(length(drugs) == length(allDrugs)) && sum(length(selectedCellLines) == length(allCellLines'))
            print(sprintf('%s DRC plot',(drugs{drug})),'-dpdf','-bestfit');
        else
            titles = [cellLineLabels{:}];
            print(sprintf('%s %s DRC plot',(drugs{drug}),titles),'-dpdf','-bestfit');
        end
        
        cd('/Users/sdalin/Dropbox (MIT)/Biology PhD/Matlab Scripts/384WellPlateReaderDRCAnalysis')
        
        
        clf
    end
            
    clf

end