%This function takes two structures with many fields (A and B) and combines
%the fields of both into A.

function [A] = combineFields(A,B)
    if isempty(A)
        A = [A;B];
    else
        drugsA = fieldnames(A);
        drugsB = fieldnames(B);
        
        for drug = 1:size(drugsB,1);
            fitsB = fieldnames(B.(drugsB{drug}));
            
            for fit = 1:size(fitsB,1);
                A.(drugsB{drug}).(fitsB{fit}) = B.(drugsB{drug}).(fitsB{fit});
            end
        end
    end
end