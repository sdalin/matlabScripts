%This function will take in PI values, cell# values, concentration values
%and a list of drug names from DRCs then fit the following equation to the
%data, output a list of drugs with IC50<10mM (including their IC50), and plot the dose response
%curves and fits for those drugs.
%pvalByConc has drugs in rows and PI values in columns, corresponding to
%concbyDrug where drugs are again in rows and the corresonding drug
%concentration is in each column.

%INPUT:
%bigstructNormed has each plate as a field, and within each field, has
%PI-values and cell#'s of the 384 well plate as a vertical vector.
%concentrations is a vertical vector of the concentrations used in each row
%of the plate.  This only works if each row had the same concentration
%across all plates.

%Equation: 
%Y=Ainf+(A0-Ainf)/(1+(X/IC50)^n)
%Ainf is expected PI-% at very high conc (ie, 0)
%A0 is expected PI-% at 0mM (ie, usually 100)
%Y is PI-%
%X is concentration
%IC50 is self explanatory
%n is Hill coefficient - measure of steepness of curve

function [fittedStruct,individualFittedStruct] = hillFitv2(bigstructNormed,concentrations)
%First define the equation we are fitting to
hill = fittype( 'Ainf+((A0-Ainf)/(1+(x/IC50)^n))', 'independent', 'x', 'dependent', 'y' );

%Change the intial values and bounds of the coeffients to help the fit
opts = fitoptions(hill);
opts.Display = 'Off';
opts.Lower = [0 -1 0 -1];
opts.StartPoint = [1 0 0.1 1];
opts.Upper = [2 0.3 1000 Inf];

%Run the actual fit for the normalized cell counts in bigstructNormed
%create structure where each plate has its own field in which columns are
%different cell lines (as on the original plate) and rows are A0,Ainf,IC50,n,rsquare,adjrsquare
plates = fieldnames(bigstructNormed);

fittedStruct = struct;
fittedHill = cell(length(bigstructNormed.(char(plates(1,:))))/16,1);
gof = cell(length(bigstructNormed.(char(plates(1,:))))/16,1);

individualFittedHill = cell(length((bigstructNormed.(char(plates(1,:))))/16)*3,1);
individualGof = cell(length((bigstructNormed.(char(plates(1,:))))/16)*3,1);

for k = 1:size(plates,1)
    fittedStruct.(char(plates(k,:))) = zeros(9,20);

    %This next for loop does the fit for each column of the 384-well plate
    for cellLine = 1:size(bigstructNormed.(char(plates(k,:))),1)/16
        %Check if there are any NANs in the normalized cell counts we're
        %about to fit and if there are, remove that value as well as the
        %concentration that goes along with it;
        thisCellLine = bigstructNormed.(char(plates(k,:)))(1+(16*(cellLine-1)):(16+(16*(cellLine-1))),1:3);
        %Make it vertical to work with the fitting script
        thisCellLineVert = reshape(thisCellLine,numel(thisCellLine),1);
        %Tile concentrations next to it
        numReps = size(thisCellLineVert,1)/16;
        currentConcs = repmat(concentrations,numReps,1);
        %Remove rows where there is a NAN in the data
        currentConcs(isnan(thisCellLineVert)) = [];
        thisCellLineVert(isnan(thisCellLineVert)) = [];
        
        %Now do the fit on the values with NANs removed:

        if length(thisCellLine) > 3
            [fittedHill{cellLine}, gof{cellLine}] = fit(currentConcs,thisCellLineVert, hill, opts );
            
            %Calculate the standard error of each variable from fitting all
            %the data
            
            %First we need the area under 95% of the t distribution given
            %reported degrees of freedom:
            alphaup = 1-0.05/2;
            upp = tinv(alphaup,gof{cellLine}.dfe);
            
            %variable to calculate sdErr for
            stErrVariables = fieldnames(fittedHill{cellLine});
            
            for field = 1:size(stErrVariables,1)
                %Now based on the reported confidence interval, we can
                %calculate the stderr (confirmed by crosschecking with fitnlm
                %calculated sterrs)
                upperConf = confint(fittedHill{cellLine});
                estimatedStandardError = (upperConf(2,field) - fittedHill{cellLine}.(char(stErrVariables(field,:))))/upp;

                %Add the sterrs to a subfield of gof called stErr.A0,
                %stErr.Ainf, etc.
                gof{cellLine}.stErr.(char(stErrVariables(field,:))) = estimatedStandardError;
            end
            
            %Do the fit on each replicate individually:
            for replicate = 1:numReps
                %Remove rows where there is a NAN in the data
                currentConcs = concentrations;
                currentCellLine = thisCellLine(:,replicate);
                currentConcs(isnan(currentCellLine(:))) = [];
                currentCellLine(isnan(currentCellLine(:))) = [];
                
                if length(currentCellLine) > 3
                    %The fit
                    [individualFittedHill{((cellLine-1)*3)+replicate}, individualGof{((cellLine-1)*3)+replicate}] = fit(currentConcs,currentCellLine, hill, opts);
                else
                    continue
                end
            end
            
            %calculate the standard error of each variable from doing three
            %fits seperately
            
            for field = 1:size(stErrVariables,1)
                %Get standard error directly from the replicates
                %Add the sterrs to a subfield of gof called stErr.A0,
                %stErr.Ainf, etc.
                
                currentVariable = zeros(3,1);
                for replicate = 1:numReps
                    if ~isempty(individualFittedHill{((cellLine-1)*3)+replicate})
                        currentVariable(numReps,1) = individualFittedHill{((cellLine-1)*3)+replicate}.(char(stErrVariables(field,:)));
                    else
                        continue
                    end
                end
                individualGof{((cellLine-1)*3)+1}.stErr.(char(stErrVariables(field,:))) = std(currentVariable)/sqrt(numReps);
            end
                

        else
            continue
        end
        
        
        %Load desired variables (best fit, SE, and rsquare) into
        %fittedStruct and individualFittedStruct
        for field = 1:size(stErrVariables,1)
            fittedStruct.(char(plates(k,:)))(field,cellLine) = fittedHill{cellLine}.(char(stErrVariables(field,:)));          
            fittedStruct.(char(plates(k,:)))(field+4,cellLine) = gof{cellLine}.stErr.(char(stErrVariables(field,:))); 
            
            for replicate = 1:numReps
                if ~isempty(individualFittedHill{((cellLine-1)*3)+replicate})
                    individualFittedStruct.(char(plates(k,:)))(field,((cellLine-1)*3)+replicate) = individualFittedHill{((cellLine-1)*3)+replicate}.(char(stErrVariables(field,:)));
                else
                    continue
                end
                
                if isfield(individualGof{((cellLine-1)*3)+replicate},'stErr')
                    individualFittedStruct.(char(plates(k,:)))(field+4,((cellLine-1)*3)+replicate) = individualGof{((cellLine-1)*3)+replicate}.stErr.(char(stErrVariables(field,:)));
                else
                    continue
                end
            end
        end
        
        fittedStruct.(char(plates(k,:)))(9,cellLine) = gof{cellLine}.rsquare;
        
        for replicate = 1:numReps
            if isfield(individualGof{((cellLine-1)*3)+replicate},'rsquare')
                individualFittedStruct.(char(plates(k,:)))(9,((cellLine-1)*3)+replicate) = individualGof{((cellLine-1)*3)+replicate}.rsquare;
            else
                continue
            end
        end
    end
    
%     %Are IC50's of each cell line and it's clones significantly different from
%     %each other? (ANOVA)
%     
%     %First need to generate random data because matlab's ANOVA requires
%     %data, not descriptive statistics
%     
%     presumedN = 43;
%     randData = zeros(presumedN,size(bigstructNormed.(char(plates(k,:))),1)/16);
%     for cellLine = 1:size(bigstructNormed.(char(plates(k,:))),1)/16
%         randData(:,cellLine) = normrnd(fittedStruct.(char(plates(k,:)))(3,cellLine),(fittedStruct.(char(plates(k,:)))(7,cellLine))*sqrt(presumedN),presumedN,1);
%     end
%     
%     %Now run the ANOVA
%     [p,tbl,stats] = anova1(randData);
    
end




%Make list of cell lines and their IC50s, stErr, and rsquare values for each plate
%then save it

rowNames = {'A0';'Ainf';'IC50';'n';'A0 SE';'Ainf SE';'IC50 SE';'n SE';'rsquare'};
cd('/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/SSC Heterogeneity DRCs/matlabOutput')
for k = 1:size(plates,1)
    cellLineIC50cell = cell(9,length(fittedStruct.(char(plates(k,:)))));
    individualCellLineIC50cell = cell(9,length(individualFittedStruct.(char(plates(k,:)))));
    
    for cellLine = 1:length(fittedStruct.(char(plates(k,:))))
        for dataVariable = 1:size(fittedStruct.(char(plates(k,:))),1)
            cellLineIC50cell{dataVariable,cellLine} = fittedStruct.(char(plates(k,:)))(dataVariable,cellLine);
        end
    end
        
    for cellLine = 1:length(individualFittedStruct.(char(plates(k,:))))
        for dataVariable = 1:size(individualFittedStruct.(char(plates(k,:))),1)
            individualCellLineIC50cell{dataVariable,cellLine} = individualFittedStruct.(char(plates(k,:)))(dataVariable,cellLine);
        end
    end
        
    T = cell2table(cellLineIC50cell);
    indT = cell2table(individualCellLineIC50cell);
    writetable(T,sprintf('%s',char(plates(k,:)),'.csv'));
    writetable(indT,sprintf('%s','individual',char(plates(k,:)),'.csv'));
end



%Plot it
% 
% for cellLine = 1:length(fieldnames(bigstructNormed))
%     hillPlot140217(reallyKillingDrugs,fittedHill,concbyDrug,pvalByConc,plate,betterHits);
% 
%     %Plot all the drugs
%     alldrugs = find(drugIC50mat(:));
%     hillPlot140217(drugnames,fittedHill,concbyDrug,pvalByConc,plate,alldrugs);

cd('/Users/sdalin/Dropbox (MIT)/Biology PhD/Matlab Scripts')





