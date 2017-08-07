thing2 = thing
thing2(find(cellfun(@isempty,thing))) = {'nan'}

for row = 1:size(thing2,1)
    for col = 1:size(thing2,2)
        thing3(row,col) = thing2{row,col};
    end
end