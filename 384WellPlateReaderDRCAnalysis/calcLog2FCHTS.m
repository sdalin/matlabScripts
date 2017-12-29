%This function will take in datAfterFit, containing fit objects and
%calculated EC50s for all the cell lines calculated thus far.  It will
%calculate log2(FC) over the parental cell line for each new experiment and
%store the data back in dataAfterFit

%INPUT:
%dataAfterFit is a structure with a field for the fit objects (subfields
%are drug names, then cell line names, then both a matrix with
%concentrations in first column, dates of measurements in first row,
%and normalized variabilities in other fields, as well as two more fields:
%individual fits with the subfields as dates of measurement, and within each
%the appropriate fit object, and lastly, a field for the fit with all the data combined,
%containing the appropriate fit object) and a field
%for paramaters of the fits (subfields are fold-changes (subfields are drug
%names, then a cell array with cell line names in first row and measured
%fold-changes below that) and EC50s (subfields are drug names, then a cell
%array with cell line names in the first row, and measured EC50s below
%that).  This file is updated with new data each time this script is run.


%OUTPUT:
%Within dataAfterFit.fitParams, within each drug's field, within the Log2FC
%field, new Log2FC values are added to the appropriate row/column of the
%struct array.  Also the 'dataLog2Calcd' field is updated with the name of each
%cellLine/Drug/Experiment as it is run.

function [dataAfterFit] = calcLog2FCHTS(dataAfterFit)

    %Value to assign to uncalculatable FC's
    uncalcFC = 20;

    %Now that all the EC50s for all the drugs and all the cell lines in
    %the new experiments have been calculated/extracted, we can calculate log2(FC) for each
    %cell line/drug over the parental.  Calculate and store in
    %appropriate section of dataAfterFit

    %First set up lists of drugs, experiments to calculate Log2FC for, and
    %cell lines
    drugs = fieldnames(dataAfterFit.rawData);
    allExperiments = dataAfterFit.fitParams.(drugs{1}).EC50(2:end,1);
    if ~isfield(dataAfterFit.fitParams.(drugs{1}),'Log2FC')
        completedExperiments = cell(0);
    elseif ~isempty(dataAfterFit.fitParams.(drugs{1}).Log2FC)
        completedExperiments = dataAfterFit.fitParams.(drugs{1}).Log2FC(2:end,1);
    else
        whut
    end
    %experiments = setdiff(allExperiments,completedExperiments);
    experiments = allExperiments;
    cellLines = fieldnames(dataAfterFit.rawData.(drugs{1}));

    %First set up log of what data has been analyzed.
    if ~isfield(dataAfterFit.fitParams,'dataLog2Calcd')
        dataAfterFit.fitParams.dataLog2Calcd = cell(0);
    end

    for experiment = 1:size(experiments,1)
        for drug = 1:size(drugs,1)
            %On first iteration, add a row to FC fields with
            %the experiment name
            if ~isfield(dataAfterFit.fitParams.(drugs{drug}),'Log2FC') %the drug has never had a log2FC calculated before
                dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC = {[];experiments{experiment}};
            elseif isempty(dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC) %A log2FC field exists but is empty
                dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC = {[];experiments{experiment}};
            elseif isfield(dataAfterFit.fitParams,sprintf('%s',drugs{drug})) && sum(strcmp(dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC(:,1),experiments{experiment}))==0 %the field exists and other log2FCs have been calculated, but not for this experiment.
                dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC = [dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC;cell(1,size(dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC,2))];
                dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC(end,1) = experiments(experiment);
            end


            columnParental = find(strcmp(dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).EC50(1,:),'Parental'));
            rowParental = find(strcmp(dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).EC50(:,1),experiments{experiment}));

            %only continue if the parental cell line has data for this
            %drug.  Otherwise note if parental is in allKiller or
            %nonKiller.
            parentalEC50 = dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).EC50(rowParental,columnParental);
            
            if isempty(cell2mat(parentalEC50))

                %Find if this drug for the parental line is in allKiller or
                %nonKiller and label it accordingly
                if sum(strcmp(dataAfterFit.fitParams.allKiller,sprintf('%s_%s_%s',drugs{drug},'Parental',experiments{experiment})))
                    parentalEC50 = {'allKiller'};
                elseif sum(strcmp(dataAfterFit.fitParams.nonKiller,sprintf('%s_%s_%s',drugs{drug},'Parental',experiments{experiment})))
                    parentalEC50 = {'nonKiller'};
                end
                
            end

            for cellLine = 1:size(cellLines,1)

                %Check if this experiment/cellLine/drug has already
                %been run
                if ~sum(ismember(dataAfterFit.fitParams.dataLog2Calcd,sprintf('%s_%s_%s',drugs{drug},cellLines{cellLine},experiments{experiment})))
                    
                    %Find if this drug for the current cell line is
                    %allKiller or nonKiller or a number and label accordingly.
                    if sum(strcmp(dataAfterFit.fitParams.allKiller,sprintf('%s_%s_%s',drugs{drug},cellLines{cellLine},experiments{experiment})))
                        currentEC50 = {'allKiller'};
                    elseif sum(strcmp(dataAfterFit.fitParams.nonKiller,sprintf('%s_%s_%s',drugs{drug},cellLines{cellLine},experiments{experiment})))
                        currentEC50 = {'nonKiller'};
                    else
                        columnCurrentEC50 = find(strcmp(dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).EC50(1,:),cellLines{cellLine}));
                        currentEC50 = dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).EC50(rowParental,columnCurrentEC50);
                    end

                    parentalEC50Mat = cell2mat(parentalEC50);
                    currentEC50Mat = cell2mat(currentEC50);
                    parentalEC50Type = whos('parentalEC50Mat');
                    currentEC50Type = whos('currentEC50Mat');
                    
                    if strcmp(parentalEC50Type.class,'char')
                        %Assign the log2FC to be 10 or -10, depending on
                        %where the parental and current EC50s are
                        if strcmp(parentalEC50,'nonKiller')
                            if strcmp(currentEC50Type.class,'double') %Parental never kills, current has an EC50 -> CS
                                Log2FC = -uncalcFC; 
                            elseif strcmp(currentEC50,'nonKiller') %Both parental and current never kill -> no change
                                Log2FC = 0;
                            elseif strcmp(currentEC50,'allKiller') %Parental never kills, current always kills -> CS
                                Log2FC = -uncalcFC;
                            end
                        elseif strcmp(parentalEC50,'allKiller')
                            if strcmp(currentEC50Type.class,'double') %Parental always kills, current has an EC50 -> CR
                                Log2FC = uncalcFC; 
                            elseif strcmp(currentEC50,'nonKiller') %Parental always kills and current never kills -> CR
                                Log2FC = uncalcFC;
                            elseif strcmp(currentEC50,'allKiller') %Parental always kills, current always kills -> no change
                                Log2FC = 0;
                            end  
                        end
                           
                    elseif strcmp(parentalEC50Type.class,'double')
                        if strcmp(currentEC50Type.class,'double') %Both parental and current have EC50s and can calculate numeric FC
                            Log2FC = log2(cell2mat(currentEC50)/cell2mat(parentalEC50));
                        elseif strcmp(currentEC50,'nonKiller') %Parental has an EC50 and current never kills -> CR
                            Log2FC = 10;
                        elseif strcmp(currentEC50,'allKiller') %Parental has an EC50, current always kills -> CS
                            Log2FC = -10;
                        end
                    end
                        
                    %Now store it
                    [rowLog2FC,columnLog2FC] = find(strcmp(dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC,sprintf('%s',cellLines{cellLine})));
                    rowExperiment = find(strcmp(dataAfterFit.fitParams.(drugs{drug}).Log2FC,experiments{experiment}));
                    %If that cell line has no data yet, add new column with
                    %cell line name and log2(FC)
                    if isempty(rowLog2FC)
                        cellLineAndLog2FC = cell(size(dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC,1),1);
                        cellLineAndLog2FC(1) = {sprintf(cellLines{cellLine})};
                        cellLineAndLog2FC(rowExperiment) = {Log2FC};
                        dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC(:,size(dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC,2)+1) = cellLineAndLog2FC;
                        %If there's already data for that cell line, add the
                        %new log2(FC) to the bottom of the column.  Also add a row
                        %to all columns, to keep dimensions acceptable for
                        %matlab
                    else
                        dataAfterFit.fitParams.(sprintf('%s',drugs{drug})).Log2FC(rowExperiment,columnLog2FC) = {Log2FC};
                    end
                end
                dataAfterFit.fitParams.dataLog2Calcd{end+1} = sprintf('%s_%s_%s',drugs{drug},cellLines{cellLine},experiments{experiment});
            end
        end
    end
end
