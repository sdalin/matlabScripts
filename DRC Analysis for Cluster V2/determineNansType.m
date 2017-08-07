[row,col] = find(isnan(heatmapMatrix2));
nanLines = cellLines(row);
nanDrugs = drugs(col);

nanList = cell(size(nanLines,1),1);
for nanCondition = 1:size(nanLines,1)
    nanList(nanCondition,1) = {sprintf('drug_%s_%s',nanDrugs{nanCondition},nanLines{nanCondition})};
end

nanList = regexprep(nanList,' ','_');

nansFit = cell(0);
nansAllKiller = cell(0);
nansNonKiller = cell(0);
nansBadFit = cell(0);

for condition = 1:size(nanList,1)
    fitted = ~isempty(find(~cellfun(@isempty,strfind(dataAfterFit.fitParams.dataFit,nanList{condition}))));
    allKiller = ~isempty(find(~cellfun(@isempty,strfind(dataAfterFit.fitParams.allKiller,nanList{condition}))));
    nonKiller = ~isempty(find(~cellfun(@isempty,strfind(dataAfterFit.fitParams.nonKiller,nanList{condition}))));
    badFit = ~isempty(find(~cellfun(@isempty,strfind(dataAfterFit.fitParams.badFit,nanList{condition}))));
    
    %check = sum(~isempty(fitted),~isempty(allKiller),~isempty(nonKiller),~isempty(badFit));
    
    fittedLogical = {[],'fitted'};
    allKillerLogical = {[],'allKiller'};
    nonKillerLogical = {[],'nonKiller'};
    badFitLogical = {[],'badFit'};
    
    %nanList(condition,2) = {sprintf('%s',badFitLogical{badFit + 1})};

    
    nanList(condition,2) = {sprintf('%s,%s,%s,%s',fittedLogical{fitted + 1},allKillerLogical{allKiller + 1},nonKillerLogical{nonKiller + 1},badFitLogical{badFit + 1})};
end
    
    