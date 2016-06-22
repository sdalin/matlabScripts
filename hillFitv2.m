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

function [fittedStruct,fittedHill,gof] = hillFitv2(bigstructNormed,concentrations)
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
fittedHill = cell(length(bigstructNormed.(char(plates(1,:))))/16,4);
gof = cell(length(bigstructNormed.(char(plates(1,:))))/16,2);
nonLinFittedHill = cell(length(bigstructNormed.(char(plates(1,:))))/16,1);

for k = 1:size(plates,1)
    fittedStruct.(char(plates(k,:))) = zeros(16,20);

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
            
            %Calculate the standard error of each variable (modify to also
            %do other variables too!) 
            alphaup = 1-0.05/2;
            upp = tinv(alphaup,gof{cellLine}.dfe);
            upperConf = confint(fittedHill{cellLine});
            estimatedStandardError = (upperConf(2,3) - fittedHill{cellLine}.IC50)/upp;
        
        else
            continue
        end
        
        fitOutput = reshape(table2array(nonLinFittedHill{cellLine}.Coefficients),numel(nonLinFittedHill{cellLine}.Coefficients),1);
        %Load desired variables into fittedStruct
        fittedStruct.(char(plates(k,:)))(:,cellLine) = fitOutput;
            
    end
end


%Make list of cell lines and their IC50s and rsquare values for each plate
%then save it

rowNames = {'A0';'Ainf';'IC50';'n';'A0 SE';'Ainf SE';'IC50 SE';'n SE';'A0 SE';'Ainf tStat';'IC50 tStat';'n tStat';'A0 pValue';'Ainf pValue';'IC50 pValue';'n pValue'};
cd('/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/SSC Heterogeneity DRCs/matlabOutput')
for k = 1:size(plates,1)
    cellLineIC50cell = cell(16,length(fittedStruct.(char(plates(k,:)))));
    for cellLine = 1:length(fittedStruct.(char(plates(k,:))))/16
        for dataVariable = 1:size(fittedStruct.(char(plates(k,:))),1)
            cellLineIC50cell{dataVariable,cellLine} = fittedStruct.(char(plates(k,:)))(dataVariable,cellLine);
%         cellLineIC50cell{2,cellLine} = fittedStruct.(char(plates(k,:)))(2,cellLine);
%         cellLineIC50cell{3,cellLine} = fittedStruct.(char(plates(k,:)))(3,cellLine);
%         cellLineIC50cell{4,cellLine} = fittedStruct.(char(plates(k,:)))(4,cellLine);
%         cellLineIC50cell{5,cellLine} = fittedStruct.(char(plates(k,:)))(5,cellLine);
%         cellLineIC50cell{6,cellLine} = fittedStruct.(char(plates(k,:)))(6,cellLine);
        end
    end
    cellLineIC50cell = [rowNames,cellLineIC50cell];
    T = cell2table(cellLineIC50cell);
    writetable(T,sprintf('%s',char(plates(k,:)),'.txt'));
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





