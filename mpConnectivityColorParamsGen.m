function [p task] = mpConnectivityColorParamsGen

KbName('UnifyKeyNames')

p.saveFile = 1;
p.testingLocation = 'laptop'; 

switch p.testingLocation
   case 'laptop' 
        p.screenSize = [12 16]; % (in)
        p.screenRes = [768 1024]; % 60 Hz. Must set screen to this!
        p.viewDist = 19.5; % (in)
        p.triggerCode = KbName('5%');
        p.keyCodes = [KbName('1!')];
    otherwise
        error('Testing location not found!')
end

% addition by AM: include the TR and make sure block lengths align with it.
p.TR = 2.25; % for 3T using previous RD ep2d_neuro_... sequence

p.condNames = {'leftM','leftP','rightM','rightP'};%,'blank'};
%p.blocksToInclude = repmat([1 1 2 2 3],1,3);
p.blocksToInclude = 1:length(p.condNames);
%p.blocksToInclude = [1 2]% 3];
p.condOrder = generateBlockSequenceColor(p.blocksToInclude);
p.numCycles = length(p.condOrder); % a cycle means one stimulus block
% This is the block length and should be an integer multiple of the TR!
p.cycleDuration = 54; % (s) - updated 2020-08-04 to increase "block" length
assert(mod(p.cycleDuration, p.TR)==0);
p.nWedgePhases = 1; % 1 wedge phase in cycle (always full screen)

p.fixSize = 0.1; 
p.orientationOption = 'set'; % 'rand','set'
switch p.orientationOption
    case 'set'
        p.gratingOrientations = 0:180/6:180-180/6;
end
p.gratingSize = []; % defaults to screen size when empty

p.maskType = 'rectangleBlur';
switch p.maskType
    case 'transparentWedge'
        p.mappingType = 'hemifieldLR';
        p.wedgeRadial = [0 30]; % radius (internal and external) of wedge (deg) eg. [0.5 12]
        p.wedgePolar = 365*(pi/180); % (rad)
    case 'transparentGaussian'
        p.maskSigmaProp = 0.1; % proportion of grating size
        p.maskAmp = 1;
    case 'rectangleBlur'
        % rectangle will be a proportion of the size of the grating stimulus
        p.apertureProp = .9;
        p.maskSmoothSize = 50;
    case 'none'
        % no special settings
    otherwise
        error('Mask type not found.')
end

p.temporalWaveform = 'sin'; % triangular, sin
% How long to stay blank before and after the stimuli (s)
% Again, these must be a multiple of the TR!
p.blankDuration = [0 0]; % 9 seconds = 4TRs of blank at the beginning
% This aligns with the paper, in which 4TRs at the beginning were discarded
assert(all(mod(p.blankDuration, p.TR)==0));
% This is added to the cycleDuration and should be an integer multiple of the TR!
%p.responseDuration = 2.25; % (s)
%assert(mod(p.responseDuration, p.TR)==0);
%p.trialCushion = 0.5; % (s) actual response window = responseDuration - trialCushion
% Finally, calculate the total length of stimulus presentation
% and verify it's a multiple of the TR
%p.total_length = (p.numCycles * (p.cycleDuration+p.responseDuration) + sum(p.blankDuration));
p.total_length = (p.numCycles * (p.cycleDuration) + sum(p.blankDuration));
assert(mod(p.total_length, p.TR)==0);

task.maxTargets = round(p.cycleDuration/4);
task.responseOptions = 1;
task.targetCushion = 1; % time cushion on each side of targets
task.ISISuggestedSigma = 2;
task.targetDur = 0.3;
task.targetAmp = 1;
task.radialBounds = [1 15]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust for each subject [M P blank]
task.unitSigma = [.2 .4 .50];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

