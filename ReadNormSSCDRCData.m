%function designed to take iQue signature data in a specified format (SSC data, each column a different clone, each plate one drug and one parental cell ine)
%one has to specificy how many rows the CSV file has.

%OUTPUT:
%bigstructNormed has normalized cell counts of each plate in alpha order.  Columns of each plate are different SSCs of the same cell line. 

function [bigstructNormed] = ReadNormSSCDRCData(filename,endRow)
    filename = '/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/SSC Heterogeneity DRCs/All CSV Files/160511_SSC_DRC_Fixed.csv';
    endRow = 1933;
    %import first two columns of spreadsheet as a cell (these have
    %information on what each well is)
    startRow = 2;
    firsttwocolumns = importfirsttwotocell(filename, startRow, endRow);
    
    %import last two columns of csv as a matrix (these are PI-% and cell #
    %respectively)
    lasttwo = importlasttwotomat(filename, startRow, endRow);
    
    %finds locations of new plates
    indexnewplate = strfind(firsttwocolumns(:,1),'Plate');
    indexnewplate = find(~cellfun(@isempty,indexnewplate));
    
    %establish plate boundaries
    platenum = length(indexnewplate);
    platebounds = zeros(platenum,2);
    ii = 1;
    while ii < platenum 
        platebounds(ii,1) = indexnewplate(ii)+2;
        platebounds(ii,2) = indexnewplate(ii+1)-1;
        ii = ii+1;
    end
    platebounds(platenum,1) = indexnewplate(platenum)+2;
    platebounds(platenum,2) = length(firsttwocolumns);

    %create structure where each plate has its own field consisting of last
    %two columns
    %first have to remove "" spaces and :
    nospaces = strrep(firsttwocolumns,' ','');
    noquotesnospaces = strrep(nospaces,'"','');
    allgone = strrep(noquotesnospaces,':','');
    k = 2;
    bigstruct = struct(allgone{indexnewplate(1)},lasttwo(platebounds(1,1):platebounds(1,2),1:2));
    while k < platenum + 1
        bigstruct = setfield(bigstruct,allgone{indexnewplate(k),1},lasttwo(platebounds(k,1):platebounds(k,2),1:2));
        k = k+1;
    end
    
    %For each plate, normalize cell # in each column to the average of rows
    %H & P, which had DMSO.
    
    %Normalize each column of each plate to average of row H & P which had
    %DMSO
    
    for k = 1:platenum
        for cellLine = 1:20
            DMSOs = [bigstruct.(allgone{indexnewplate(k),1})(8 + (16 * (cellLine - 1)),2),bigstruct.(allgone{indexnewplate(k),1})(16 + (16 * (cellLine - 1)),2)];
            DMSOmean = mean(DMSOs);
            bigstructNormed.(allgone{indexnewplate(k),1})(1 + (16 * (cellLine - 1)):16 + (16 * (cellLine - 1)),1) = bigstruct.(allgone{indexnewplate(k),1})(1 + (16 * (cellLine - 1)):16 + (16 * (cellLine - 1)),2)./DMSOmean;
        end
    end
        

           
% 
%     %Extract all cell numbers and concentrations.  Each row of these two
%     %matrices is equivalent to an entire plate.
%     k = 2;
%     while k < platenum + 1
%         cellNum = bigstruct.(allgone{indexnewplate(k),1})(:,2);
%         MLSconc = bigstruct.(allgone{indexnewplate(k),1})(:,1);
%         k = k + 1
%     end
% 
%    
%    
%     %Extract drugnames with more than cellcut cells - for DRCs we want all
%     %wells
%     cellcut=0;
%     LDlog = cellNum > cellcut;
%     
%     namelog = logical([0;0;0;LDlog]);
%     cellLineNames = allgone(namelog,3); %cant use allgone for names because it has header
%     cellLineNames = unique(cellLineNames);
%     
%     plate = allgone{indexnewplate(1)};
%     plate2 = plate;
%     plate = plate(20:end);
%     if isempty(plate)
%         plate=plate2;
%     end
%     %Find out what concentrations were used on this plate (MODIFIED FOR REF
%     %ONLY)
%     %concentrations = zeros(input(sprintf('Input #concentrations tested on plate %s', plate))); 
%     concentrations = zeros(8);
%     
%     %Find the PI values and cell number and concentrations of each drug in drugnames for each of the
%     %concentrations specified
%     pvalByConc = zeros(length(cellLineNames),length(concentrations));
%     cellNumByConc = zeros(length(cellLineNames),length(concentrations));
%     concbyDrug = zeros(length(cellLineNames),length(concentrations));
%     for drug = 1:length(cellLineNames)
%         %first get the index of where this drug is
%         drugIndex = strcmp(allgone(indexMLSplate+2:platebounds(find(platebounds(:,1)==indexMLSplate+2),2),3),cellLineNames(drug));
%         %Next find the concentrations this drug was tested with
%         drugConcs = sort(MLSconc(drugIndex),'ascend');
%         concbyDrug(drug,:) = drugConcs';
%         %Now find the row of this drug and each concentration
%         for conc = 1:length(drugConcs)
%             concentrationIndex = bigstruct.(allgone{indexMLSplate(1),1})(:,1) == drugConcs(conc);
%             index = find(concentrationIndex & drugIndex);
%             pvalByConc(drug,conc) = bigstruct.(allgone{indexMLSplate,1})(index,2);
%             cellNumByConc(drug,conc) = bigstruct.(allgone{indexMLSplate,1})(index,4);
%         end
%     end
    