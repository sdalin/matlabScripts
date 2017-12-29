%function designed to take Tecan plate reader DRC data along with a 'drugs'
%tab which indicates the drug in each well, and a 'dose' tab that indicates
%the dose in each well.  Also requires reference file to translate drug
%codes to actual drug names

%INPUT:filename is the name of the file that contains data.  Sheets with
%raw data must be named 'Sheet#' and sheet with concentration/drug name
%info must be named 'Dose' and 'Drugs'.  Cell line on each plate must be in
%E12.  Should be a new file for each replicate. No drug controls must be
%columns 2 and 24.

%OUTPUT:
%bigstructNormed has normalized viability of each plate in a field named by 
%the cell line name, subfield 'rawData'.  Other subfield 'CV' is the CV of
%the no drug wells. 
%Info field has drug name and dose for every well of the plate.
%Columns of each plate are different drugs on the same cell line. 
%Normalization  divides by median of no-media wells (rezasurin salt data)

function [bigstructNormed] = ReadNormPlateReaderHTSDRCData(filename,cellLinePlateCompoundBarcodes,Directory)

%% Start by extracting the meta-data and raw-data from the excel files
    %Find number of sheets in the excel file 'filename'
    [status,sheets] = xlsfinfo(filename);
    numOfSheets = numel(sheets);
    
    %for each experiment, import raw fluorescence data and cell line name from all data sheets of filename .xls file
    %and import the drug and concentration information into a field
    %called 'Info'
    for sheet = 1:numOfSheets
        if strcmp(sheets{sheet},'Dose')
            %Grab the doses
            [num,text,doses] = xlsread(filename,'Dose','B24:Y39');
            
            %Store doses in bigstructNormed.Info
            bigstructNormed.Info.doses = doses;
            
        elseif strcmp(sheets{sheet},'Drugs')
            %Grab the drugNames
            [num,drugNames,raw] = xlsread(filename,'Drugs','B24:Y39');
            
            %Store drug names in BigstructNormed.Info
            bigstructNormed.Info.drugNames = drugNames;
            
        else
            %Grab the raw data from this sheet
            [num,text,cellRawFluorescence] = xlsread(filename, sprintf('%s',sheets{sheet}), 'B2:17');
            removeOVERs = cellRawFluorescence;
            removeOVERs(find(strcmp(cellRawFluorescence,'OVER'))) = {NaN};    
            rawFluorescence = cell2mat(removeOVERs);
            %Grab the cell line name from this sheet
            [num,cellLine,raw] = xlsread(filename,sheet,'E12');
            cellLine = strrep(cellLine,' ','_');
            cellLine = strrep(cellLine,',','');
            bigstruct.(sprintf('%s',char(cellLine))) = rawFluorescence;
        end
    end

%% Remove background fluorescence (no need with luminescence?  Or at least, there are no wells without cells!    
%     %For each plate, remove background fluorescence in each column by
%     %subtraction the median of column 23, which had no cells, store in bigstructNormed.
%     for cellLine = 1:size(fieldnames(bigstruct),1)
%         cellLines = fieldnames(bigstruct);
%         backgroundFluorescence = nanmedian(bigstruct.(sprintf('%s',cellLines{cellLine}))(:,23));
%         bigstructNoBackground.(sprintf('%s',cellLines{cellLine})) = bigstruct.(sprintf('%s',cellLines{cellLine}))-backgroundFluorescence;
%     end

%% Normalize viability to median of DMSO-only wells (columns 2 & 24)

    %Normalize each column of each plate to average of columns 2 and 24
    %which had no drugs
    cellLines = fieldnames(bigstruct);
    for cellLine = 1:size(fieldnames(bigstruct),1)
        noDrugs = [bigstruct.(sprintf('%s',cellLines{cellLine}))(:,2),bigstruct.(sprintf('%s',cellLines{cellLine}))(:,24)];
        medianNoDrug = nanmedian(noDrugs(:));
        CV = std(noDrugs(:))/mean(noDrugs(:));
        bigstructNormed.(sprintf('%s',cellLines{cellLine})).rawData = bigstruct.(sprintf('%s',cellLines{cellLine}))/medianNoDrug;
        bigstructNormed.(sprintf('%s',cellLines{cellLine})).CV = CV;
    end
    
end
   
