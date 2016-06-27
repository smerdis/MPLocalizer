function [seq count] = generateBlockSequenceColor(blocks)

% 1 is M, 2 is P, 3 is blank
% blocks = [1 1 1 2 2 2 3 3 3];

count = 0;
seqOK = 0;
while ~seqOK
    order = randperm(length(blocks));
    seq = blocks(order);
    seqOK = 1; % until found otherwise
    
    % test whether the same condition follows itself
    if any(diff(seq)==0)
        seqOK = 0; % failed the repeat test
    end
    
    count = count+1;
end
