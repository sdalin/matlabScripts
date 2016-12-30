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

%Equation:
%Y=Ainf+(A0-Ainf)/(1+(X/IC50)^n)
%Ainf is expected PI-% at very high conc (ie, 0)
%A0 is expected PI-% at 0mM (ie, usually 100)
%Y is PI-%
%X is concentration
%IC50 is self explanatory
%n is Hill coefficient - measure of steepness of curve

function [dataAfterFit] = hillFitv2(bigstructNormed,dataAfterFit,folderName)
    %First define the equation we are fitting to
    hill = fittype( 'Ainf+((A0-Ainf)/(1+(x/IC50)^n))', 'independent', 'x', 'dependent', 'y' );

    %Change the intial values and bounds of the coeffients to help the fit
    opts = fitoptions(hill);
    opts.Display = 'Off';
    opts.Lower = [0.7 -0.3 0 0];
    opts.StartPoint = [1 0 1 1];
    opts.Upper = [1.5 0.3 10000 10];

    %First set up log of what data has been analyzed.
    if ~isfield(dataAfterFit.fitParams,'dataFit')
        dataAfterFit.fitParams.dataFit = cell(0);
    end


    %Run the actual fit for the normalized viabilities in bigstructNormed
    experiments = fieldnames(bigstructNormed);
    for experiment = 1:size(experiments,1)
        cellLines = fieldnames(bigstructNormed.(experiments{experiment}));
        for cellLine = 1:size(cellLines,1)
            if ~strcmp(cellLines{cellLine},'Info')
                drugs = bigstructNormed.(experiments{experiment}).Info.drugNames;
                for drug = 1:size(drugs,2)
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

                    if isempty(drugs{drug})
                        continue

                        %Check if this experiment/cellLine/drug has already
                        %been run.  Don't enter into analysis if it already exists

                    elseif dataNotFit
                        %Check if there are any NANs in the normalized viability we're
                        %about to fit and if there are, remove that value as well as the
                        %concentration that goes along with it;
                        thisCellLine = bigstructNormed.(sprintf('%s',experiments{experiment})).(sprintf('%s',cellLines{cellLine}))(:,drug+2);
                        %Extract current concentrations to a vector
                        currentConcs = bigstructNormed.(sprintf('%s',experiments{experiment})).Info.doses(:,drug);
                        %Remove rows where there is a NAN in the data
                        currentConcs(isnan(thisCellLine)) = [];
                        thisCellLine(isnan(thisCellLine)) = [];

                        %Now do the fit on the values with NANs removed, and re-do fit if its bad (user-input):

                        fitGood = 'no';
                        while strcmp(fitGood,'no')
                            %check if there is already a graph of this cell
                            %line in the data folder, and if yes, move on
                            %to the next drug.
                            if length(thisCellLine) > 3
                                [fittedHill{cellLine}, gof{cellLine}] = fit(currentConcs,thisCellLine,hill,opts);

                                %Run a script to plot the fit and check if it looks
                                %good.  Then another script to accept user input on
                                %whether to keep or throw away all/some of the
                                %data.

                                %Make the plot
                                plotIndividualCellLines(fittedHill{cellLine},currentConcs,thisCellLine,drugs{drug},cellLines{cellLine},experiments{experiment},folderName);

                                %Have user input if the data is good, bad, or
                                %what data to remove
                                fitGood = input('Is the fit and data good? (enter "yes" or "no")','s');
                                if strcmp(fitGood,'yes')
                                    %save the fit object and gof into dataAfterFit
                                    dataAfterFit.fitObj.IndFit.(sprintf('drug_%s',drugs{drug})).(sprintf('%s_%s',cellLines{cellLine},experiments{experiment})) = fittedHill{cellLine};
                                    dataAfterFit.gofObj.IndFit.(sprintf('drug_%s',drugs{drug})).(sprintf('%s_%s',cellLines{cellLine},experiments{experiment})) = gof{cellLine};
                                    skipToNextCellLine = 0;



                                    %Save the data
                                    save('dataAfterFit','dataAfterFit');

                                    %Then continue to next iteration of while
                                    %loop, which will exit while loop and move
                                    %on in the code
                                    continue
                                else
                                    removeData = input('What data to remove? (enter "all" or concentration values in a vector)','s');
                                    if strcmp(removeData,'all')
                                        %Delete the plot, and do not save any
                                        %of the data, then continue
                                        delete(sprintf('%s/matlabOutput/%s %s %s.pdf',folderName,drugs{drug},cellLines{cellLine},experiments{experiment}));
                                        fittedHill{cellLine} = [];
                                        gof{cellLine} = [];
                                        skipToNextCellLine = 1;

                                        %Save the data
                                        save('dataAfterFit','dataAfterFit');

                                        break
                                    else
                                        %Delete the specified data points from currentConcs and thisCellLine,
                                        %remove the plot
                                        removeData = str2num(removeData);
                                        thisCellLine(find(ismember(currentConcs,removeData))) = [];
                                        currentConcs(find(ismember(currentConcs,removeData))) = [];
                                        delete(sprintf('%s/matlabOutput/%s %s %s.pdf',folderName,drugs{drug},cellLines{cellLine},experiments{experiment}));                                    skipToNextCellLine = 0;
                                        skipToNextCellLine = 0;

                                        %Save the data
                                        save('dataAfterFit','dataAfterFit');
                                    end
                                end

                            else
                                %If there are less than 3 datapoints, cannot do
                                %the fit.  Just move on with the next drug/cellLine.
                                break
                            end
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
                            cellLineAndEC50(1:2) = {sprintf(cellLines{cellLine});currentEC50};
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
                    end
                end
            end
        end
    end
end




