function [p, t, task] = runMPConnectivityBehavPilot(subjectID)

condNames = {'leftM','leftP','rightM','rightP'};
% make the block sequence subject to constraints
[seq, ~] = generateBlockSequenceColor([1 1 1 2 2 2 3 3 3 4 4 4]);
for i=1:length(seq)
    [p, t, task] = runMPConnectivity(subjectID, seq(i), condNames{seq(i)}) ;
end