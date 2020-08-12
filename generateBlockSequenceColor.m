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
    elseif length(unique(blocks)) > 3 &&  (any(diff(seq)==2) || any(diff(seq)==-2))
        % we have >3 conditions, meaning connectivity (leftM, leftP, rightM, rightP)
        % we don't want leftM to be followed by rightM, etc
        seqOK = 0;
    end
    
    count = count+1;
end
