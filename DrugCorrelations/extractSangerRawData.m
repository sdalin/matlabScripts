%Extract raw dose and viability data from Sanger data

function [rawData] = extractSangerRawData(directory)

    %Import the data
    dataFilename = 'v17a_public_raw_data.csv';
    last92 = repmat('%s',1,92);
    formatSpec = ['%d%s%{dd-MMM-yyyy}D%d%d%s%d%s%d%s',last92];
    rawImport = readtable(sprintf('%s%s',directory,dataFilename));
    rawData = [rawImport(:,1),rawImport(:,9:end)];
    
    %For each drug and cell line, extract median max viability, median
    %control viability, normalized treated viability (subtract control,
    %divide by max), and the doses it was tested against.
    
    %Import drug and cell line ID to names
    drugFilename = 'Screened_Compounds.xlsx';
    drugInfo = readtable(sprintf('%s%s',directory,drugFilename));
   
%% this should actually loop over IC_results_ID to take into account any replicates...
    
    normedData = rawData(:,1:15);
    normedData.Properties.VariableNames(7:end) = [{'normed_max'},{'normed2'},{'normed3'},{'normed4'},{'normed5'},{'normed6'},{'normed7'},{'normed8'},{'normed9'}];
    
    for ICresult = 1:size(rawData,1)
        drugName = drugInfo.DrugName(drugInfo.DrugID==rawData.DRUG_ID(ICresult));
        
        %Make the drug name an acceptable field name
        %Fix any drug names so it can be a fieldname
        if isnumeric(drugName) %For drug names that are entirely numbers
            drugName = strcat('drug_',int2str(drugName));
        else
            drugName = ['drug_',drugName];
        end
        drugName = strrep(drugName,'-','_');
        drugName = strrep(drugName,' ','_');
        drugName = strrep(drugName,'(','');
        drugName = strrep(drugName,')','');
        drugName = strrep(drugName,'[','');
        drugName = strrep(drugName,']','');
        drugName = strrep(drugName,'/','_');
        drugName = strrep(drugName,'\','_');
                
        cellLine = rawData.CELL_LINE_NAME(ICresult);
        
        %Extract control and blank fluorescence
        colNames = rawData.Properties.VariableNames;
        controlCols = cellfun(@(x) ~isempty(x),regexp(colNames,'control\d+'));  
        controlDouble = str2double(rawData{ICresult,controlCols});
        control = nanmedian(controlDouble);
        
        blankCols = cellfun(@(x) ~isempty(x),regexp(colNames,'blank\d+'));
       	blankDouble = str2double(rawData{ICresult,blankCols});
        blank = nanmedian(blankDouble);
        
        %Now extract and normalize the treated fluorescent data
        rawCols = cellfun(@(x) ~isempty(x),regexp(colNames,'raw\d+|raw_max'));
        rawDouble = str2double(rawData{ICresult,rawCols});
        
        rawNoBackground = rawDouble - blank;
        rawNormalized = rawNoBackground/control;
        
        %Save into normedData table under raw
        normedData(ICresult,7:end) = num2cell(rawNormalized);
        
        
        
    end
        
end


%readTable is reading some numeric columns as chars.  Setting a formatspec
%for the import doesn't work for some reason, so instead just detect the
%issues and fix them here:
function [fixedArray] = correctDataType(badArray,ICresult,indices)
    fixedArray = zeros(size(indices));
    for index = 1:size(indices,2)
        if ~isnumeric(badArray{ICresult,indices(index)})
            fixedArray(1,index) = str2double(badArray{ICresult,indices(index)});
        else
            fixedArray(1,index) = badArray{ICresult,indices(index)};
        end
    end
end
