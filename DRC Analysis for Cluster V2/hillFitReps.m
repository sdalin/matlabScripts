function [dataAfterFit] = hillFitReps(dataAfterFit)
    %First define the equation we are fitting to
    hill = fittype( 'Ainf+((A0-Ainf)/(1+(x/IC50)^n))', 'independent', 'x', 'dependent', 'y' );

    %Change the intial values and bounds of the coeffients to help the fit
    opts = fitoptions(hill);
    opts.Display = 'Off';
    opts.Lower = [0.7 -0.3 0 0];
    opts.StartPoint = [1 0 1 1];
    opts.Upper = [1.5 0.3 10000 10];
    opts.Exclude = [];    

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

end