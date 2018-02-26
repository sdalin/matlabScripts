%Wrapper script for fitting Sanger raw data into 4-parameter hill fit,
%extracting EC50s and AUCs.

%First Extract the data as a table and save in a .mat file
directory = '/Users/sdalin/Dropbox (MIT)/Biology PhD/2018/H&L Labs/DrugSensitivityCorrelations/Data/';

[rawData] = extractSangerRawData(directory);

%Now do the fit to 4 parameter hill curve.  save in a .mat file.

%Plot 10 random cell lines per drug

%Extract EC50s and AUCs from the hill fit.