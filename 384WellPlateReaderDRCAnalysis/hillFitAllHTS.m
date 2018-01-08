%This function will take in normalized viability, concentration values,
%and drug names from DRCs in a big structure, then fits the following equation to the
%data, outputs fit objects, raw EC50s, and EC50 FC's in a new struct.

%INPUT:
%bigstructNormed has each day of data collection from the current folder as a field.
%Within those, each plate is a field, and there is also an 'Info' field containing drug
%names and concentrations of drug in each well on the plates in that sub-field.
%Plate fields have normalized viability in 384-well plate format.
%dataAfterFit is a saved file containing all data previously generated,
%which will be updated and added to by running this script.

%OUTPUT:
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
%different cell line names. This file is updated with new data each time this script is run.
%Also makes a cell array will the names of all drugs in the current folder
%in the first row and all cell lines in the current folder in the second
%row.  Saves the variable to the matlabOutput folder as
%drugsCellLinesThisFolder

%Equation:
%Y=Ainf+(A0-Ainf)/(1+(X/IC50)^n)
%Ainf is expected PI-% at very high conc (ie, 0)
%A0 is expected PI-% at 0mM (ie, usually 100)
%Y is PI-%
%X is concentration
%IC50 is self explanatory
%n is Hill coefficient - measure of steepness of curve

function [drugsCellLinesThisFolder] = hillFitAllHTS(bigstructNormed,folderName)
    %First define the equation we are fitting to
    hill = fittype( 'Ainf+((A0-Ainf)/(1+(x/IC50)^n))', 'independent', 'x', 'dependent', 'y' );

    %Change the intial values and bounds of the coeffients to help the fit
    opts = fitoptions(hill);
    opts.Display = 'Off';
    opts.Lower = [0.7 -0.3 0 0];
    opts.StartPoint = [1 0 1 1];
    opts.Upper = [1.5 0.3 10000 10];
    opts.Exclude = [];

    %Initialize dataAfterFit
    dataAfterFit = struct('rawData',[],'fitObj',[],'gofObj',[],'fitParams',[]);
    dataAfterFit.fitObj.IndFit = [];
    dataAfterFit.fitObj.allDataFit = [];
    dataAfterFit.gofObj.IndFit = [];
    dataAfterFit.gofObj.allDataFit = [];
    
    %Set up log of what data has been analyzed.
    if ~isfield(dataAfterFit.fitParams,'dataFit')
        dataAfterFit.fitParams.dataFit = cell(0);
    end

    if ~isfield(dataAfterFit.fitParams,'allKiller')
        dataAfterFit.fitParams.allKiller = cell(0);
    end

    if ~isfield(dataAfterFit.fitParams,'nonKiller')
        dataAfterFit.fitParams.nonKiller = cell(0);
    end
    
    if ~isfield(dataAfterFit.fitParams,'lessThanThreeDataPoints')
        dataAfterFit.fitParams.lessThanThreeDataPoints = cell(0);
    end  
    
    if ~isfield(dataAfterFit.fitParams,'badFit')
        dataAfterFit.fitParams.badFit = cell(0);
    end
    
    if ~isfield(dataAfterFit.fitParams,'allDataFit')
        dataAfterFit.fitParams.allDataFit = cell(0);
    end
    
    drugsCellLinesThisFolder = cell(2,0);
    
    %Initialize EC50 fields
    experiments = fieldnames(bigstructNormed);
    for experiment = 1:size(experiments,1)
        if strcmp(experiments{experiment},'Info')
            continue
        elseif strcmp(experiments{experiment},'Directory')
            continue
        end
        
        drugs = unique(bigstructNormed.(experiments{experiment}).Info.drugNames);
        
        %remove 'DMSO' and 'EMPTY' from drug names
        drugs = drugs(~strcmp(drugs,'DMSO'));
        drugs = drugs(~strcmp(drugs,'EMPTY'));
        
        
        %Initialize EC50 fields with
        %the experiment names in first column and cell lines in
        %first row
        for drug = 1:size(drugs,1)
            
            %Find actual drug name
            if sum(strcmp(drugs{drug},bigstructNormed.Info(:,1)))
                wheredrug = strcmp(drugs{drug},bigstructNormed.Info(:,1));
                %Remove any charachters that make the drug name an invalid
                %field name
                drugName = bigstructNormed.Info(wheredrug,2);
                drugName = strrep(drugName,' ','_');
                drugName = strrep(drugName,'-','');
                drugName = strrep(drugName,',','');
                drugName = strrep(drugName,'(','');
                drugName = strrep(drugName,')','');
                drugName = char(drugName);
            else
                drugName = drugs{drug};     
            end
            
            
            if isfield(dataAfterFit.fitParams,sprintf('drug_%s',drugName)) && sum(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugName)).EC50(:,1),experiments{experiment}))==0
                dataAfterFit.fitParams.(sprintf('drug_%s',drugName)).EC50 = [dataAfterFit.fitParams.(sprintf('drug_%s',drugName)).EC50;cell(1,size(dataAfterFit.fitParams.(sprintf('drug_%s',drugName)).EC50,2))];
                dataAfterFit.fitParams.(sprintf('drug_%s',drugName)).EC50(end,1) = experiments(experiment);
            elseif ~isfield(dataAfterFit.fitParams,sprintf('drug_%s',drugName))
                dataAfterFit.fitParams.(sprintf('drug_%s',drugName)).EC50 = {[];experiments{experiment}};
            end

        end
    end
    
    save('dataAfterFit','dataAfterFit');
    clear dataAfterFit
    
    %Run the actual fit for the normalized viabilities in bigstructNormed
    experiments = fieldnames(bigstructNormed);
    for experiment = 1:size(experiments,1)
        if strcmp(experiments{experiment},'Info')
            continue
        elseif strcmp(experiments{experiment},'Directory')
            continue
        end
        
        cellLines = fieldnames(bigstructNormed.(experiments{experiment}));
        
        %Initialize some variables for temporary data storage
        allKiller = cell(0);
        nonKiller = cell(0);
        badFit = cell(0);
        dataFit = cell(0);
        lessThanThreeDataPoints = cell(0);
        fitObj = struct;
        gofObj = struct;
        rawData = struct;
        fitParams = struct; 
        
        for cellLine = 1:size(cellLines,1)
            if ~(strcmp(cellLines{cellLine},'Info')) && ~(strcmp(cellLines{cellLine},'Directory'))
                
                %Indicate this cell line is in this folder
                if sum(strcmp(drugsCellLinesThisFolder(2,:),cellLines{cellLine})) < 1
                    drugsCellLinesThisFolder{2,end+1} = cellLines{cellLine};
                end
                
                %List all drugs tested on this cell line
                drugs = unique(bigstructNormed.(experiments{experiment}).Info.drugNames);
                
                %remove 'DMSO' and 'EMPTY' from drug names
                drugs = drugs(~strcmp(drugs,'DMSO'));
                drugs = drugs(~strcmp(drugs,'EMPTY'));                
                
                for drug = 1:size(drugs,1)
                   
                    experiment
                    drug
                    cellLine
                    
                    %Find actual drug name
                    if sum(strcmp(drugs{drug},bigstructNormed.Info(:,1)))
                        wheredrug = strcmp(drugs{drug},bigstructNormed.Info(:,1));
                        %Remove any charachters that make the drug name an invalid
                        %field name
                        drugName = bigstructNormed.Info(wheredrug,2);
                        drugName = strrep(drugName,' ','_');
                        drugName = strrep(drugName,'-','');
                        drugName = strrep(drugName,',','');
                        drugName = strrep(drugName,'(','');
                        drugName = strrep(drugName,')','');
                        drugName = char(drugName);
                    else
                        drugName = drugs{drug};
                    end
                    
                    %Find cells that this drug is in
                    currentDrugLocation = strcmp(bigstructNormed.(sprintf('%s',experiments{experiment})).Info.drugNames,drugs{drug});
                    
                    %Indicate this drug is in this folder
                    if sum(strcmp(drugsCellLinesThisFolder(1,:),drugs{drug})) < 1
                        drugsCellLinesThisFolder{1,end+1} = drugs{drug};
                    end
                                                    
                    if isempty(drugs{drug})
                        continue
                    end
               
                    %Check if there are any NANs in the normalized viability we're
                    %about to fit and if there are, remove that value as well as the
                    %concentration that goes along with it;
                    thisCellLine = bigstructNormed.(sprintf('%s',experiments{experiment})).(cellLines{cellLine}).rawData(currentDrugLocation);
                    %Extract current concentrations to a vector
                    currentConcs = bigstructNormed.(sprintf('%s',experiments{experiment})).Info.doses(currentDrugLocation);
                    currentConcs = cell2mat(currentConcs);
                    %Remove rows where there is a NAN in the data
                    currentConcs(isnan(thisCellLine)) = [];
                    thisCellLine(isnan(thisCellLine)) = [];

                    %Now do the fit on the values with NANs removed, and re-do fit if its bad (residuals 3sd away from mean):


                    %Check if all viability values are >0.7 (nonKiller) or
                    %<0.3 (allKiller).  Indicate which of those
                    %it falls into, save the cell line/drug name into 
                    %the allKiller or nonKiller fitParams field and do not 
                    %save any of the data, then continue

                    fittedHill{cellLine} = [];
                    fittedHillNoOutliers{cellLine} = [];
                    gof{cellLine} = [];
                    outliers = [];

                    %Check if at least one point has viability greater
                    %than 0.2
                    if sum(thisCellLine > 0.2) < 1
                        skipToNextDrug = 1;
                        allKiller{end+1} = sprintf('drug_%s_%s_%s',drugName,cellLines{cellLine},experiments{experiment});

                        %make scatterplot and save to
                        %'rejectedFits' folder
                        plotIndividualCellLines(fittedHill{cellLine},fittedHillNoOutliers{cellLine},currentConcs,thisCellLine,outliers,'not fit',drugs{drug},cellLines{cellLine},experiments{experiment},folderName,skipToNextDrug);

                    %Check if at least one point has viability less
                    %than 0.8
                    elseif sum(thisCellLine < 0.8) < 1
                        skipToNextDrug = 1;
                        nonKiller{end+1} = sprintf('drug_%s_%s_%s',drugName,cellLines{cellLine},experiments{experiment});

                        %make scatterplot and save to
                        %'rejectedFits' folder
                        plotIndividualCellLines(fittedHill{cellLine},fittedHillNoOutliers{cellLine},currentConcs,thisCellLine,outliers,'not fit',drugs{drug},cellLines{cellLine},experiments{experiment},folderName,skipToNextDrug);

                    elseif length(thisCellLine) > 3
                        %clear previous cell line's excluded data
                        opts.Exclude = [];

                        %fit the data to the model
                        [fittedHill{cellLine}, gof{cellLine}] = fit(currentConcs,thisCellLine,hill,opts);

                        %Now check for residuals of fit and remove
                        %any >3 sd from the fit, remove those
                        %points and re-make the fit.  
                        fdata = feval(fittedHill{cellLine},currentConcs);
                        residuals = fdata - thisCellLine;
                        outliers = abs(residuals) > 3*std(residuals); 
                        outliersIndicies = excludedata(currentConcs,thisCellLine,'indices',outliers);

                        opts.Exclude = outliersIndicies;
                        [fittedHillNoOutliers{cellLine}, gofNoOutliers{cellLine}] = fit(currentConcs,thisCellLine,hill,opts);

                            
%                         %Check that Rsquare is greater than threshold
%                         %of 0.8 and exclude if not
%                         if gofNoOutliers{cellLine}.rsquare < 0.8
%                             skipToNextDrug = 1;
% 
%                             badFit{end+1} = sprintf('drug_%s_%s_%s',drugName,cellLines{cellLine},experiments{experiment});
% 
%                             %make scatterplot and save to
%                             %'rejectedFits' folder
%                             plotIndividualCellLines(fittedHill{cellLine},fittedHillNoOutliers{cellLine},currentConcs,thisCellLine,outliers,gofNoOutliers{cellLine}.rsquare,drugName,cellLines{cellLine},experiments{experiment},folderName,skipToNextDrug); 
%                         end

                        skipToNextDrug = 0;

                        %Make and save the plot of this cell line
                        plotIndividualCellLines(fittedHill{cellLine},fittedHillNoOutliers{cellLine},currentConcs,thisCellLine,outliers,gofNoOutliers{cellLine}.rsquare,drugName,cellLines{cellLine},experiments{experiment},folderName,skipToNextDrug);

                        %save the fit object and gof into dataAfterFit
                        fitObj.IndFit.(sprintf('drug_%s',drugName)).(sprintf('%s_%s',cellLines{cellLine},experiments{experiment})) = fittedHillNoOutliers{cellLine};
                        gofObj.IndFit.(sprintf('drug_%s',drugName)).(sprintf('%s_%s',cellLines{cellLine},experiments{experiment})) = gofNoOutliers{cellLine};

                        %Indicate this drug/cell line was fit
                        dataFit{end+1} = sprintf('drug_%s_%s_%s',drugName,cellLines{cellLine},experiments{experiment});

                        %indicate data checking has happened
                        checkData = true;

                    else
                        %If there are less than 3 datapoints, cannot do
                        %the fit.  Just move on with the next drug/cellLine.
                        skipToNextDrug = 1;
                        lessThanThreeDataPoints{end+1} = sprintf('drug_%s_%s_%s',drugName,cellLines{cellLine},experiments{experiment});

                        %make scatterplot and save to
                        %'rejectedFits' folder
                        plotIndividualCellLines(fittedHill{cellLine},fittedHillNoOutliers{cellLine},currentConcs,thisCellLine,outliers,'not fit',drugName,cellLines{cellLine},experiments{experiment},folderName,skipToNextDrug);

                        %break

                        clf
                    end


                    %Now save the raw data into another field of this mega
                    %struct.  First check if there's any raw data there yet and
                    %if not start the field with this raw data and if yes, add the new data
                    %vertically below whats already in there. %Need to nest
                    %the if statements to be able to check if the drug
                    %field has the cell line subfield.  (errors if there is
                    %no drug field, but what we really want is inclusive
                    %or...)

                    rawData.(sprintf('drug_%s',drugName)).(cellLines{cellLine}) = [num2cell(currentConcs),num2cell(thisCellLine)];
                    
                    %Now extract the EC50 of each cell line/drug into the fitParams field of dataAfterFit
                    if ~skipToNextDrug
                        currentEC50 = fitObj.IndFit.(sprintf('drug_%s',drugName)).(sprintf('%s_%s',cellLines{cellLine},experiments{experiment})).IC50;
                   
                    else
                        currentEC50 = [];
                    end
                    
                    %Store the EC50 into the row/column corresponding to the
                    %experiment/cellLine
                    
                    if ~isfield(fitParams,sprintf('drug_%s',drugName))
                        fitParams.(sprintf('drug_%s',drugName)).EC50 = cell(0);
                    end
                    
                    [row,column] = find(strcmp(fitParams.(sprintf('drug_%s',drugName)).EC50,sprintf('%s',cellLines{cellLine})));
                    
                    %If that cell line has no data yet, add new column with
                    %cell line name and EC50
                    if isempty(row)
                        cellLineAndEC50 = cell(size(fitParams.(sprintf('drug_%s',drugName)).EC50,1),1);
                        row = 2;
                        cellLineAndEC50(1) = {sprintf(cellLines{cellLine})};
                        cellLineAndEC50(row) = {currentEC50};
                        fitParams.(sprintf('drug_%s',drugName)).EC50(:,size(fitParams.(sprintf('drug_%s',drugName)).EC50,2)+1) = cellLineAndEC50;
                    end
                    
                end
            end
        end
        
        java.lang.System.gc()
        
        %Now need to load all this data back into dataAfterFit
        load('dataAfterFit','dataAfterFit');
        
        %Load notations of which cell line/drug pairs are which quality
        dataAfterFit.fitParams.allKiller = [dataAfterFit.fitParams.allKiller,allKiller];
        dataAfterFit.fitParams.nonKiller = [dataAfterFit.fitParams.nonKiller,nonKiller];
        dataAfterFit.fitParams.badFit = [dataAfterFit.fitParams.badFit,badFit];
        dataAfterFit.fitParams.dataFit = [dataAfterFit.fitParams.dataFit,dataFit];
        dataAfterFit.fitParams.lessThanThreeDataPoints = [dataAfterFit.fitParams.lessThanThreeDataPoints,lessThanThreeDataPoints];
        
        %Load fit data
        dataAfterFit.fitObj.IndFit = combineFields(dataAfterFit.fitObj.IndFit,fitObj.IndFit);    
        dataAfterFit.gofObj.IndFit = combineFields(dataAfterFit.gofObj.IndFit,gofObj.IndFit);
        
        %Load raw data
        drugs = fieldnames(rawData);
        for drug = 1:size(drugs,1)
            cellLinesToAdd = fieldnames(rawData.(drugs{drug}));
            
            if isfield(dataAfterFit.rawData,drugs{drug})
                cellLinesInDataAfterFit = fieldnames(dataAfterFit.rawData.(drugs{drug}));
            else
                cellLinesInDataAfterFit={};
            end
            
            for cellLine = 1:size(cellLinesToAdd,1)
                
                if strcmp(cellLines{cellLine},'Info')
                    continue
                end
                
                if isempty(cellLinesInDataAfterFit) && ~isfield(rawData.(drugs{drug}),cellLines{cellLine})
                    continue
                elseif isempty(cellLinesInDataAfterFit) && isfield(rawData.(drugs{drug}),cellLines{cellLine})
                    dataAfterFit.rawData.(drugs{drug}).(sprintf('%s',cellLines{cellLine})) = rawData.(drugs{drug}).(cellLines{cellLine});
                elseif ~sum(strcmp(cellLines{cellLine},cellLinesInDataAfterFit)) && ~isfield(rawData.(drugs{drug}),cellLines{cellLine})
                    continue
                elseif ~sum(strcmp(cellLines{cellLine},cellLinesInDataAfterFit)) && isfield(rawData.(drugs{drug}),cellLines{cellLine})
                    dataAfterFit.rawData.(drugs{drug}).(sprintf('%s',cellLines{cellLine})) = rawData.(drugs{drug}).(cellLines{cellLine});
                elseif isfield(dataAfterFit.rawData,drugs{drug}) && ~isfield(rawData.(drugs{drug}),cellLines{cellLine})
                    continue
                elseif isfield(dataAfterFit.rawData,drugs{drug}) && isfield(rawData.(drugs{drug}),cellLines{cellLine})
                    dataAfterFit.rawData.(drugs{drug}).(sprintf('%s',cellLines{cellLine})) = ([dataAfterFit.rawData.(drugs{drug}).(sprintf('%s',cellLines{cellLine}));rawData.(drugs{drug}).(cellLines{cellLine})]);
                else
                    pause
                end
            end
            
        end
        
        %Load EC50 data
        drugs = fieldnames(fitParams);
        
        for drug = 1:size(drugs,1)
            cellLines = fitParams.(drugs{drug}).EC50(1,1:end);
            
            for cellLine = 1:size(cellLines,2)
                [row,column] = find(strcmp(dataAfterFit.fitParams.(drugs{drug}).EC50,(cellLines{cellLine})));
                
                if isempty(row)
                    [row,tmp] = find(strcmp(dataAfterFit.fitParams.(drugs{drug}).EC50,(experiments{experiment})));
                    dataAfterFit.fitParams.(drugs{drug}).EC50(1,end+1) = {sprintf(cellLines{cellLine})};
                    dataAfterFit.fitParams.(drugs{drug}).EC50(row,end) = fitParams.(drugs{drug}).EC50(2,cellLine);
                else
                    [row,tmp] = find(strcmp(dataAfterFit.fitParams.(drugs{drug}).EC50,(experiments{experiment})));
                    dataAfterFit.fitParams.(drugs{drug}).EC50(row,column) = fitParams.(drugs{drug}).EC50(2,cellLine);
                end
            end
        end
        
        %Save the data
        save('dataAfterFit','dataAfterFit');
        
        %Clear that variable
        clear dataAfterFit
        
    end     

    load('dataAfterFit');
    %Now that we have all the data ever collected for this cell
    %line all together for this cell line, do the hill fit on ALL
    %the data
    
    drugs = fieldnames(dataAfterFit.rawData);
    for drug = 1:size(drugs,1)
        cellLines = fieldnames(dataAfterFit.rawData.(drugs{drug}));
        for cellLine = 1:size(cellLines,1)
            if ~strcmp(cellLines{cellLine},'Info')
            
                if ~isfield(dataAfterFit.rawData.(drugs{drug}),(sprintf('%s',cellLines{cellLine})))
                    continue
                end
    
                opts.Exclude = [];
                [fittedHill{cellLine}, gof{cellLine}] = fit(cell2mat(dataAfterFit.rawData.(drugs{drug}).(sprintf('%s',cellLines{cellLine}))(:,1)),cell2mat(dataAfterFit.rawData.(drugs{drug}).(sprintf('%s',cellLines{cellLine}))(:,2)),hill,opts);
   
                %and save the fit object and gof into dataAfterFit
                dataAfterFit.fitObj.allDataFit.(drugs{drug}).(sprintf('%s',cellLines{cellLine})) = fittedHill{cellLine};
                dataAfterFit.gofObj.allDataFit.(drugs{drug}).(sprintf('%s',cellLines{cellLine})) = gof{cellLine};
            end
        end
    end

    %save dataAfterFit
    save('dataAfterFit','dataAfterFit');
    clear dataAfterFit
    
    %save drugs and cell lines in this folder
    cd(sprintf('%s/matlabOutput',folderName))
    save('drugsCellLinesThisFolder','drugsCellLinesThisFolder')
    cd(sprintf('%s../',folderName))

end






