
%1 runs all data, 2 through 5 run dox, vin, pac, cis, in that order
for i = 1:5
    dataType = {'AllData','Dox','Vin','Pac','Cis'};
    [eva{i},differenceReconData{i},differenceCoeffsPCAPPCA{i},differenceCoeffsALSPPCA{i},differenceCoeffsALSPCA{i},evaPCA{i},evaPPCA{i},evaALS{i}] = Log2EC50Analysis(dataAfterFit,folder,dataType{i});
end
    
