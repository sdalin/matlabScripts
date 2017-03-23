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

function [dataAfterFit,drugsCellLinesThisFolder] = hillFitv3(bigstructNormed,dataAfterFit,folderName)
    tic;
    %First define the equation we are fitting to
    hill = fittype( 'Ainf+((A0-Ainf)/(1+(x/IC50)^n))', 'independent', 'x', 'dependent', 'y' );

    %Change the intial values and bounds of the coeffients to help the fit
    opts = fitoptions(hill);
    opts.Display = 'Off';
    opts.Lower = [0.7 -0.3 0 0];
    opts.StartPoint = [1 0 1 1];
    opts.Upper = [1.5 0.3 10000 10];
    opts.Exclude = [];

    %First set up log of what data has been analyzed.
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
    
    drugsCellLinesThisFolder = cell(2,0);
    
    %Run the actual fit for the normalized viabilities in bigstructNormed
    experiments = fieldnames(bigstructNormed);
    for experiment = 1:size(experiments,1)
        cellLines = fieldnames(bigstructNormed.(experiments{experiment}));
        for cellLine = 1:size(cellLines,1)
            if ~strcmp(cellLines{cellLine},'Info')
                
                %Indicate this cell line is in this folder
                if sum(strcmp(drugsCellLinesThisFolder(2,:),cellLines{cellLine})) < 1
                    drugsCellLinesThisFolder{2,end+1} = cellLines{cellLine};
                end
                
                drugs = bigstructNormed.(experiments{experiment}).Info.drugNames;
                for drug = 1:size(drugs,2)
                    
                    %Indicate this drug is in this folder
                    if sum(strcmp(drugsCellLinesThisFolder(1,:),drugs{drug})) < 1
                        drugsCellLinesThisFolder{1,end+1} = drugs{drug};
                    end
                    
                    %On first iteration, add a row to EC50 fields with
                    %the experiment name
                    if isfield(dataAfterFit.fitParams,sprintf('drug_%s',drugs{drug})) && sum(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(:,1),experiments{experiment}))==0
                        dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50 = [dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50;cell(1,size(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50,2))];
                        dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(end,1) = experiments(experiment);
                    elseif ~isfield(dataAfterFit.fitParams,sprintf('drug_%s',drugs{drug}))
                        dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50 = {[];experiments{experiment}};
                    end

                    %This variable is only true if the data has never before
                    %been fit & saved
                    dataNotFit = ~sum(strcmp(dataAfterFit.fitParams.dataFit,sprintf('drug_%s_%s_%s',drugs{drug},cellLines{cellLine},experiments{experiment})));
                    dataNotAllKiller = ~sum(strcmp(dataAfterFit.fitParams.allKiller,sprintf('drug_%s_%s_%s',drugs{drug},cellLines{cellLine},experiments{experiment})));
                    dataNotNonKiller = ~sum(strcmp(dataAfterFit.fitParams.nonKiller,sprintf('drug_%s_%s_%s',drugs{drug},cellLines{cellLine},experiments{experiment})));
                    dataNotBad = dataNotAllKiller && dataNotNonKiller;
                                           
                    if isempty(drugs{drug})
                        continue

                        %Check if this experiment/cellLine/drug has already
                        %been run.  Don't enter into analysis if it already exists

                    elseif dataNotFit && dataNotBad
                        %Check if there are any NANs in the normalized viability we're
                        %about to fit and if there are, remove that value as well as the
                        %concentration that goes along with it;
                        thisCellLine = bigstructNormed.(sprintf('%s',experiments{experiment})).(sprintf('%s',cellLines{cellLine}))(:,drug+2);
                        %Extract current concentrations to a vector
                        currentConcs = bigstructNormed.(sprintf('%s',experiments{experiment})).Info.doses(:,drug);
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
                        
                        %Check if at least two points are viability greater
                        %than 0.3
                        if sum(thisCellLine > 0.3) < 2
                            skipToNextCellLine = 1;
                            dataAfterFit.fitParams.allKiller{end+1} = sprintf('drug_%s_%s_%s',drugs{drug},cellLines{cellLine},experiments{experiment});

                            %make scatterplot and save to
                            %'rejectedFits' folder
                            plotIndividualCellLines(fittedHill{cellLine},fittedHillNoOutliers{cellLine},currentConcs,thisCellLine,outliers,drugs{drug},cellLines{cellLine},experiments{experiment},folderName);

                        %Check if at least two points are viability less
                        %than 0.7
                        elseif sum(thisCellLine < 0.7) < 2
                            skipToNextCellLine = 1;
                            dataAfterFit.fitParams.nonKiller{end+1} = sprintf('drug_%s_%s_%s',drugs{drug},cellLines{cellLine},experiments{experiment});

                            %make scatterplot and save to
                            %'rejectedFits' folder
                            plotIndividualCellLines(fittedHill{cellLine},fittedHillNoOutliers{cellLine},currentConcs,thisCellLine,outliers,drugs{drug},cellLines{cellLine},experiments{experiment},folderName);

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

                            %Make and save the plot of this cell line
                            plotIndividualCellLines(fittedHill{cellLine},fittedHillNoOutliers{cellLine},currentConcs,thisCellLine,outliers,drugs{drug},cellLines{cellLine},experiments{experiment},folderName);

                            %save the fit object and gof into dataAfterFit
                            dataAfterFit.fitObj.IndFit.(sprintf('drug_%s',drugs{drug})).(sprintf('%s_%s',cellLines{cellLine},experiments{experiment})) = fittedHill{cellLine};
                            dataAfterFit.gofObj.IndFit.(sprintf('drug_%s',drugs{drug})).(sprintf('%s_%s',cellLines{cellLine},experiments{experiment})) = gof{cellLine};
                            skipToNextCellLine = 0;

                            %Save the data
                            save('dataAfterFit','dataAfterFit');

                            %indicate data checking has happened
                            checkData = true;

                        else
                            %If there are less than 3 datapoints, cannot do
                            %the fit.  Just move on with the next drug/cellLine.
                            dataAfterFit.fitParams.lessThanThreeDataPoints{end+1} = sprintf('drug_%s_%s_%s',drugs{drug},cellLines{cellLine},experiments{experiment});
                            break
                            
                        clf
                        end
                    else
                        %If execution ends up in here, its because the data was
                        %already analyzed.
                        skipToNextCellLine = 1;
                    end


                    %Now save the raw data into another field of this mega
                    %struct.  First check if there's any raw data there yet and
                    %if not start the field with this raw data and if yes, add the new data
                    %vertically below whats already in there. %Need to nest
                    %the if statements to be able to check if the drug
                    %field has the cell line subfield.  (errors if there is
                    %no drug field, but what we really want is inclusive
                    %or...)
                    if ~skipToNextCellLine && dataNotFit
                        if isfield(dataAfterFit.rawData,sprintf('drug_%s',drugs{drug}))
                            if isfield(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})),sprintf('%s',cellLines{cellLine}))
                                dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})) = ([dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine}));num2cell(currentConcs),num2cell(thisCellLine)]);
                            else
                                dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})) = [currentConcs];
                                dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})) = num2cell(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})));
                                dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine}))(:,size(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})),2)+1) = ([num2cell(thisCellLine)]);
                            end
                            %This second part is repetative, but only because we want the same thing to happen if the drug isn't a field and
                            %if the cell line isn't a field, but if the drug isn't
                            %a field we can't ceck if the cell line isn't a field.
                        else
                            dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})) = [currentConcs];
                            dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})) = num2cell(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})));
                            dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine}))(:,size(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})),2)+1) = ([num2cell(thisCellLine)]);
                        end
                    end

                    %Now that we have all the data ever collected for this cell
                    %line all together for this cell, do the hill fit on ALL
                    %the data

                    if (isfield(dataAfterFit.rawData,sprintf('drug_%s',drugs{drug})) && isfield(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})),sprintf('drug_%s',drugs{drug}))) && length(cell2mat(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine}))(:,1))) < 3
                        drugs{drug}
                        cellLines{cellLine}
                        experiments{experiment}
                        cell2mat(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine}))(:,1))
                        pause;
                    end
                    
                    if ~skipToNextCellLine
                        [fittedHill{cellLine}, gof{cellLine}] = fit(cell2mat(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine}))(:,1)),cell2mat(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine}))(:,2)),hill,opts);

                        %and save the fit object and gof into dataAfterFit
                        dataAfterFit.fitObj.allDataFit.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})) = fittedHill{cellLine};
                        dataAfterFit.gofObj.allDataFit.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})) = gof{cellLine};

                        %Now extract the EC50 of each cell line/drug into the fitParams field of dataAfterFit
                        currentEC50 = dataAfterFit.fitObj.IndFit.(sprintf('drug_%s',drugs{drug})).(sprintf('%s_%s',cellLines{cellLine},experiments{experiment})).IC50;
                    else
                        currentEC50 = [];
                    end

                    %Store the EC50 into the row/column corresponding to the
                    %experiment/cellLine if this data has not been stored
                    %before
                    if dataNotFit
                        [row,column] = find(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50,sprintf('%s',cellLines{cellLine})));
                        %If that cell line has no data yet, add new column with
                        %cell line name and EC50
                        if isempty(row)
                            cellLineAndEC50 = cell(size(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50,1),1);
                            row = find(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(:,1),experiments{experiment}));
                            cellLineAndEC50(1) = {sprintf(cellLines{cellLine})};
                            cellLineAndEC50(row) = {currentEC50};
                            dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(:,size(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50,2)+1) = cellLineAndEC50;
                            %If there's already data for that cell line, add the
                            %new EC50 to the row with the current experiment
                        else
                            row = find(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(:,1),experiments{experiment}));
                            if isempty(row)
                                save('dataAfterFit','dataAfterFit');
                                'somehow the row for this experiment was not initialized'
                            else
                                dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(row,column) = {currentEC50};
                            end
                        end
                        %Indicate that this cell line/drug is done
                        dataAfterFit.fitParams.dataFit{end+1} = sprintf('drug_%s_%s_%s',drugs{drug},cellLines{cellLine},experiments{experiment});
                        save('dataAfterFit','dataAfterFit');
                    end
                end
            end
        end
    end
    
    %save drugs and cell lines in this folder
    cd(sprintf('%s/matlabOutput',folderName))
    save('drugsCellLinesThisFolder','drugsCellLinesThisFolder')
    cd('/Users/sdalin/Dropbox (MIT)/Biology PhD/Matlab Scripts/384WellPlateReaderDRCAnalysis');
    
    
end




