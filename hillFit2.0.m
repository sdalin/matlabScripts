%This function will take in PI values, cell# values, concentration values
%and a list of drug names from DRCs then fit the following equation to the
%data, output a list of drugs with IC50<10mM (including their IC50), and plot the dose response
%curves and fits for those drugs.
%pvalByConc has drugs in rows and PI values in columns, corresponding to
%concbyDrug where drugs are again in rows and the corresonding drug
%concentration is in each column.

%Equation: 
%Y=Ainf+(A0-Ainf)/(1+(X/IC50)^n)
%Ainf is expected PI-% at very high conc (ie, 0)
%A0 is expected PI-% at 0mM (ie, usually 100)
%Y is PI-%
%X is concentration
%IC50 is self explanatory
%n is Hill coefficient - measure of steepness of curve

function [outputdose,killing] = hillFit2.0(bigstructNormed)
%First define the equation we are fitting to
hill = fittype( 'Ainf+((A0-Ainf)/(1+(x/IC50)^n))', 'independent', 'x', 'dependent', 'y' );

%Change the intial values and bounds of the coeffients to help the fit
opts = fitoptions(hill);
opts.Display = 'Off';
opts.Lower = [0 -1 0 -1];
opts.StartPoint = [100 0 0.1 1];
opts.Upper = [200 100 1000 Inf];

%Run the actual fit for the normalized cell counts in bigstructNormed
%create structure where each plate has its own field in which columns are
%different cell lines (as on the original plate) and columns are A0,Ainf,IC50,n,rsquare,adjrsquare
plates = fieldnames(bigstructNormed);

fittedStruct = struct;
for k = 1:size(plates,1)
    [fittedHill{drug}, gof{drug}] = fit(concbyDrug(drug,:)',pvalByConc(drug,:)', hill, opts );
    fittedStruct = setfield(fittedStruct,plates(k),[fittedHill{drug}, gof{drug}]);
end




fittedHill = cell(length(pvalByConc),1);
gof = cell(length(pvalByConc),1);
for drug = 1:size(drugnames)
    [fittedHill{drug}, gof{drug}] = fit(concbyDrug(drug,:)',pvalByConc(drug,:)', hill, opts );
end

%Make list of drugs and their IC50s and rsquare values
drugIC50cell = cell(length(drugnames),5);
for drug = 1:size(drugnames)
    drugIC50cell{drug,1} = fittedHill{drug}.IC50;
    drugIC50cell{drug,2} = gof{drug}.rsquare;
    drugIC50cell{drug,3} = fittedHill{drug}.Ainf;
    drugIC50cell{drug,4} = fittedHill{drug}.A0;
    drugIC50cell{drug,5} = fittedHill{drug}.n;
end


%Plot it
% hillPlot140217(reallyKillingDrugs,fittedHill,concbyDrug,pvalByConc,plate,betterHits);

%Plot all the drugs
%alldrugs = find(drugIC50mat(:));
%hillPlot140217(drugnames,fittedHill,concbyDrug,pvalByConc,plate,alldrugs);


% %Add drug names to drugIC50cell for nice output
cd('/Users/sdalin/Dropbox (MIT)/Biology PhD/2015/Hemann Lab/Expts for Bos UROP/Calculating alpha/Output')
drugIC50cell = [drugnames,drugIC50cell];
T = cell2table(drugIC50cell,'VariableNames',{'DrugName','IC50','Rsquared','Ainf','A0','hillCoeff'});
writetable(T,sprintf('%s',plate,hairpin,'.txt'));





