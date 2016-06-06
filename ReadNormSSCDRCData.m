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

function [bigstruct,drugnames,pvalByConc,cellNumByConc,concbyDrug,plate] = DRCs150420(filename,endRow,hairpin)
    filename = '/Users/sdalin/Dropbox (MIT)/Biology PhD/2016/Hemann Lab/CR.CS/SSC Heterogeneity DRCs/All CSV Files/160511_SSC_DRC_Fixed.csv';
    endRow = 1933;
    %import first two columns of spreadsheet as a cell (these have
    %information on what each well is)
    startRow = 2;
    firstthreecolumns = importfirsttwotocell(filename, startRow, endRow);
    
    %import last two columns of csv as a matrix (these are PI-% and cell #
    %respectively)
    startRow = startRow + 2;
    lastsix = importlastsixtomat(filename, startRow, endRow);
    lastsix(isnan(lastsix))= 0; %so that NaN's don't mess things up
    lastsix(:,[2,3])=[]; %Remove the unneeded two columns
    
    %finds locations of new plates
    indexnewplate = strfind(firstthreecolumns(:,1),'Plate');
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
    platebounds(platenum,2) = length(firstthreecolumns);

    %create structure where each plate has its own field consisting of last
    %three columns
    %first have to remove "" spaces and :
    nospaces = strrep(firstthreecolumns,' ','');
    noquotesnospaces = strrep(nospaces,'"','');
    allgone = strrep(noquotesnospaces,':','');
    k = 2;
    bigstruct = struct(allgone{indexnewplate(1)},lastsix(platebounds(1,1):platebounds(1,2),1:4));
    while k < platenum + 1
        bigstruct = setfield(bigstruct,allgone{indexnewplate(k),1},lastsix(platebounds(k,1):platebounds(k,2),1:4));
        k = k+1;
    end
    
    %find location of hairpin plate in question.  Everywhere it says 'MLS'
    %below really refers to the hairpin in question
    
    %indexMLSplate is only the location of the plate name in the entire CSV
    indexMLSplate = strfind(firstthreecolumns(:,1),char(hairpin));
    %this next line tests if each entry of indexMLSplate (which has no
    %entries exept the row containing the word 'MLS' in firstthreecolumns)
    %is empty and puts a 1 if its got something, then find outputs that one
    %index, and that is where the MLS plate starts.
    indexMLSplate = find(~cellfun(@isempty,indexMLSplate));
    
    %finds locations of untreated wells
    %first indexes for each type of untreated the concatenates them then
    %sorts them
    indexDMSOplate = strfind(firstthreecolumns(:,2),'DMSO');
    indexDMSOplate = find(~cellfun(@isempty,indexDMSOplate));
    indexWaterplate = strfind(firstthreecolumns(:,2),'Water');
    indexWaterplate = find(~cellfun(@isempty,indexWaterplate));
    indexEtohplate = strfind(firstthreecolumns(:,2),'EtOH');
    indexEtohplate = find(~cellfun(@isempty,indexEtohplate));
    indexPosplate = strfind(firstthreecolumns(:,2),'Positive');
    indexPosplate = find(~cellfun(@isempty,indexPosplate));
    indexUT = [indexDMSOplate;indexWaterplate;indexEtohplate;indexPosplate];
    sortedindexUT = sortrows(indexUT);
    %sortedindexUT is for the entire csv file

    %Extract MLS cell numbers and concentrations
    MLScellnum = bigstruct.(allgone{indexMLSplate(1),1})(:,4);
    MLSconc = bigstruct.(allgone{indexMLSplate,1})(:,1);

    
    %identifies where in the plate are the non-UT wells
    UTlowerindMLS = (sortedindexUT <= indexMLSplate+2+length(bigstruct.(allgone{indexMLSplate(1),1})(:,2))); %cuts off higher indices
    UTplateindMLS = (sortedindexUT(UTlowerindMLS) >= indexMLSplate+2); %then cuts lower indices
    UTwellsMLSentirecsv = sortedindexUT(UTplateindMLS);
    UTwellsjustMLS = UTwellsMLSentirecsv-(indexMLSplate+1);
    nonUTlog = true(length(bigstruct.(allgone{indexMLSplate(1),1})(:,1)),1);
    nonUTlog(UTwellsjustMLS) = 0; %nonUTlog is a logical array of non UT well location
    %where UT wells are 0 (false) and treatment wells are 1 (true)
   
    %Extract drugnames with more than cellcut cells - for DRCs we want all
    %wells
    cellcut=0;
    LDlog = MLScellnum > cellcut;
    LDlog = and(LDlog,nonUTlog);
    
    namelog = logical([0;0;0;LDlog]);
    drugnames = allgone(namelog,3); %cant use allgone for names because it has header
    drugnames = unique(drugnames);
    
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
    pvalByConc = zeros(length(drugnames),length(concentrations));
    cellNumByConc = zeros(length(drugnames),length(concentrations));
    concbyDrug = zeros(length(drugnames),length(concentrations));
    for drug = 1:length(drugnames)
        %first get the index of where this drug is
        drugIndex = strcmp(allgone(indexMLSplate+2:platebounds(find(platebounds(:,1)==indexMLSplate+2),2),3),drugnames(drug));
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
    