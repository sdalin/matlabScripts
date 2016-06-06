%function designed to take iQue signature data in a specified format from CSV where
%each hairpin is represented once.  one has to specificy the LD limits, lowest cell number
%as well as how many rows the CSV file has. also make sure all hairpins are
%in caps and that the wells are all in the same order for each plate

%true? 
%need to fix the fact that for plates like S2 and S3 there are vrey
%few cells in untreated wells.  thus some wells with around 10 cells or so
%have 0 GFP+ cells which isnt taken into the weighted average when it
%should be

%should have an output that is the average of the number of cells in each
%well across each plate and consider using that average as the minimum cell
%number cut off (or maybe variance of cell number as well)

%OUTPUT:
%drugnames is array 
%bigstruct has each hairpin plate in alpha order (MLS last).  Columns of each plate are:  drug conc, PI-%, GFP+%, cell number, 
%MLS has the pivals, cell num and conc for drugs on the MLS plate.

%/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/SSC Heterogeneity DRCs/All CSV Files

function [bigstruct,cellLineNames,pvalByConc,cellNumByConc,concbyDrug,plate] = DRCs150420(filename,endRow,hairpin)
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
    
    k = 2;
    while k < platenum + 1
        %here write a for loop or a while loop that counts off each
        %384-well column and divides each value by the average.
        
        %Another possibility would be to export each plate (each leaf of
        %the structure) into a matrix, then re-shape into 16x24, then
        %normalize then re-shape again into column then store into
        %bigstructNormed.
        bigstructNormed.(allgone{indexnewplate(k),1}) = bigstruct.(allgone{indexnewplate(k),1})(:,2) / 
        
    

    %Extract all cell numbers and concentrations.  Each row of these two
    %matrices is equivalent to an entire plate.
    k = 2;
    while k < platenum + 1
        cellNum = bigstruct.(allgone{indexnewplate(k),1})(:,2);
        MLSconc = bigstruct.(allgone{indexnewplate(k),1})(:,1);
        k = k + 1
    end

   
   
    %Extract drugnames with more than cellcut cells - for DRCs we want all
    %wells
    cellcut=0;
    LDlog = cellNum > cellcut;
    
    namelog = logical([0;0;0;LDlog]);
    cellLineNames = allgone(namelog,3); %cant use allgone for names because it has header
    cellLineNames = unique(cellLineNames);
    
    plate = allgone{indexnewplate(1)};
    plate2 = plate;
    plate = plate(20:end);
    if isempty(plate)
        plate=plate2;
    end
    %Find out what concentrations were used on this plate (MODIFIED FOR REF
    %ONLY)
    %concentrations = zeros(input(sprintf('Input #concentrations tested on plate %s', plate))); 
    concentrations = zeros(8);
    
    %Find the PI values and cell number and concentrations of each drug in drugnames for each of the
    %concentrations specified
    pvalByConc = zeros(length(cellLineNames),length(concentrations));
    cellNumByConc = zeros(length(cellLineNames),length(concentrations));
    concbyDrug = zeros(length(cellLineNames),length(concentrations));
    for drug = 1:length(cellLineNames)
        %first get the index of where this drug is
        drugIndex = strcmp(allgone(indexMLSplate+2:platebounds(find(platebounds(:,1)==indexMLSplate+2),2),3),cellLineNames(drug));
        %Next find the concentrations this drug was tested with
        drugConcs = sort(MLSconc(drugIndex),'ascend');
        concbyDrug(drug,:) = drugConcs';
        %Now find the row of this drug and each concentration
        for conc = 1:length(drugConcs)
            concentrationIndex = bigstruct.(allgone{indexMLSplate(1),1})(:,1) == drugConcs(conc);
            index = find(concentrationIndex & drugIndex);
            pvalByConc(drug,conc) = bigstruct.(allgone{indexMLSplate,1})(index,2);
            cellNumByConc(drug,conc) = bigstruct.(allgone{indexMLSplate,1})(index,4);
        end
    end
    