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

%Equation:
%Y=Ainf+(A0-Ainf)/(1+(X/IC50)^n)
%Ainf is expected PI-% at very high conc (ie, 0)
%A0 is expected PI-% at 0mM (ie, usually 100)
%Y is PI-%
%X is concentration
%IC50 is self explanatory
%n is Hill coefficient - measure of steepness of curve

function [dataAfterFit] = hillFitv2(bigstructNormed,dataAfterFit)
    %First define the equation we are fitting to
    hill = fittype( 'Ainf+((A0-Ainf)/(1+(x/IC50)^n))', 'independent', 'x', 'dependent', 'y' );

    %Change the intial values and bounds of the coeffients to help the fit
    opts = fitoptions(hill);
    opts.Display = 'Off';
    opts.Lower = [0.7 -0.3 0 0];
    opts.StartPoint = [1 0 1 1];
    opts.Upper = [1.5 0.3 10000 10];

    %Run the actual fit for the normalized viabilities in bigstructNormed
    experiments = fieldnames(bigstructNormed);
    for experiment = 1:size(experiments,1)
        cellLines = fieldnames(bigstructNormed.(sprintf('%s',experiments{experiment})));
        for cellLine = 1:size(cellLines,1)
            if strcmp(cellLines{cellLine},'Info')
                continue
            else
                drugs = bigstructNormed.(sprintf('%s',experiments{experiment})).Info.drugNames;
                for drug = 1:size(drugs,2)
                    %Check if there are any NANs in the normalized viability we're
                    %about to fit and if there are, remove that value as well as the
                    %concentration that goes along with it;
                    thisCellLine = bigstructNormed.(sprintf('%s',experiments{experiment})).(sprintf('%s',cellLines{cellLine}))(:,drug+2);
                    %Extract current concentrations to a vector
                    currentConcs = bigstructNormed.(sprintf('%s',experiments{experiment})).Info.doses(:,drug);
                    %Remove rows where there is a NAN in the data
                    currentConcs(isnan(thisCellLine)) = [];
                    thisCellLine(isnan(thisCellLine)) = [];

                    %Now do the fit on the values with NANs removed:
                    if length(thisCellLine) > 3
                        [fittedHill{cellLine}, gof{cellLine}] = fit(currentConcs,thisCellLine,hill,opts);

                        %and save the fit object and gof into dataAfterFit
                        dataAfterFit.fitObj.IndFit.(sprintf('drug_%s',drugs{drug})).(sprintf('%s_%s',cellLines{cellLine},experiments{experiment})) = fittedHill{cellLine};
                        dataAfterFit.gofObj.IndFit.(sprintf('drug_%s',drugs{drug})).(sprintf('%s_%s',cellLines{cellLine},experiments{experiment})) = gof{cellLine};
                    else
                        continue
                    end

                    %Now save the raw data into another field of this mega
                    %struct.  First check if there's any raw data there yet and
                    %if not start the field with this raw data and if yes, add the new data
                    %vertically below whats already in there. %Need to nest
                    %the if statements to be able to check if the drug
                    %field has the cell line subfield.  (errors if there is
                    %no drug field, but what we really want is inclusive
                    %or...)
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
                    

                    %Now that we have all the data ever collected for this cell
                    %line all together for this cell, do the hill fit on ALL
                    %the data
                    [fittedHill{cellLine}, gof{cellLine}] = fit(cell2mat(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine}))(:,1)),cell2mat(dataAfterFit.rawData.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine}))(:,2)),hill,opts);

                    %and save the fit object and gof into dataAfterFit
                    dataAfterFit.fitObj.allDataFit.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})) = fittedHill{cellLine};
                    dataAfterFit.gofObj.allDataFit.(sprintf('drug_%s',drugs{drug})).(sprintf('%s',cellLines{cellLine})) = gof{cellLine};
                    
                    %Now extract the EC50 of each cell line/drug into the fitParams field of dataAfterFit
                    currentEC50 = dataAfterFit.fitObj.IndFit.(sprintf('drug_%s',drugs{drug})).(sprintf('%s_%s',cellLines{cellLine},experiments{experiment})).IC50;
                    %Find out if this drug has any data yet - if not, start
                    %new field.  If it does, store data in that field.
                    if isfield(dataAfterFit.fitParams,sprintf('drug_%s',drugs{drug}))
                        [row,column] = find(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50,sprintf('%s',cellLines{cellLine})));
                        %If that cell line has no data yet, add new column with
                        %cell line name and EC50
                        if isempty(row)
                            cellLineAndEC50 = cell(size(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50,1),1);
                            cellLineAndEC50(1:2) = {sprintf(cellLines{cellLine});currentEC50};
                            dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(:,size(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50,2)+1) = cellLineAndEC50;
                        %If there's already data for that cell line, add the
                        %new EC50 to the bottom of the column.  Also add a row
                        %to all columns, to keep dimensions acceptable for
                        %matlab
                        else
                            emptyCell = find(cellfun('isempty',dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(:,column)));
                            if isempty(emptyCell)
                                dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50 = [dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50;cell(1,size(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50,2))];
                                dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(size(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50,1),column) = {currentEC50};
                            else
                                dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(emptyCell,column) = {currentEC50};
                            end
                        end
                    else
                        dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50 = {sprintf(cellLines{cellLine});currentEC50};
                    end
                end
            end
        end
        
        
        
        %Now that all the EC50s for all the drugs and all the cell lines in
        %this experiment have been extracted, we can calculate log2(FC) for each
        %cell line/drug over the parental.  Calculate and store in
        %appropriate section of dataAfterFit
        for drug = 1:size(drugs,2)
            if isempty(drugs{drug})
                continue
            else
                columnParental = find(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(1,:),'Parental'));
                parentalEC50 = dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(end,columnParental);
                for cellLine = 1:size(cellLines,1)
                    if strcmp(cellLines{cellLine},'Info')
                        continue
                    else
                        columnCurrentEC50 = find(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(1,:),cellLines{cellLine}));
                        currentEC50 = dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).EC50(end,columnCurrentEC50);
                        Log2FC = log2(cell2mat(currentEC50)/cell2mat(parentalEC50));
                        %Now store it
                        %Find out if this drug has any data yet - if not, start
                        %new field.  If it does, store data in that field.
                        if isfield(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})),'Log2FC')
                            [row,column] = find(strcmp(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC,sprintf('%s',cellLines{cellLine})));
                            %If that cell line has no data yet, add new column with
                            %cell line name and log2(FC)
                            if isempty(row)
                                cellLineAndLog2FC = cell(size(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC,1),1);
                                cellLineAndLog2FC(1:2) = {sprintf(cellLines{cellLine});Log2FC};
                                dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC(:,size(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC,2)+1) = cellLineAndLog2FC;
                                %If there's already data for that cell line, add the
                                %new log2(FC) to the bottom of the column.  Also add a row
                                %to all columns, to keep dimensions acceptable for
                                %matlab
                            else
                                emptyCell = find(cellfun('isempty',dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC(:,column)));
                                if isempty(emptyCell)
                                    dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC = [dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC;cell(1,size(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC,2))];
                                    dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC(size(dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC,1),column) = {Log2FC};
                                else
                                    dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC(emptyCell,column) = {Log2FC};
                                end
                            end
                        else
                            dataAfterFit.fitParams.(sprintf('drug_%s',drugs{drug})).Log2FC = {sprintf(cellLines{cellLine});Log2FC};
                        end
                    end
                end
            end
        end 
    end
end




