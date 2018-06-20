function [p task] = mpLocalizerColorParamsGen

p.saveFile = 1;
p.testingLocation = 'laptop'; 

switch p.testingLocation
   case 'laptop' 
        p.screenSize = [12 16]; % (in)
        p.screenRes = [768 1024]; % 60 Hz. Must set screen to this!
        p.viewDist = 19.5; % (in)
        p.triggerCode = KbName('5%');
        p.keyCodes = [KbName('1!') KbName('2@') KbName('3#') KbName('4$')];
    otherwise
        error('Testing location not found!')
end

p.condNames = {'M','P','blank'};
p.blocksToInclude = repmat([1 1 2 2 3],1,3);
% p.blocksToInclude = [1 2 3];
p.condOrder = generateBlockSequenceColor(p.blocksToInclude);
p.numCycles = length(p.condOrder); % a cycle means one stimulus block
p.cycleDuration = 18; % (s)
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
p.blankDuration = [8 8];
p.responseDuration = 2; % (s)
p.trialCushion = 0.5; % (s) actual response window = responseDuration - trialCushion

task.maxTargets = 3; % maxTargets+1 is the number of intervals
task.responseOptions = 0:task.maxTargets;
task.targetCushion = 1; % time cushion on each side of targets
task.ISISuggestedSigma = 2;
task.targetDur = 0.3;
task.targetAmp = 1;
task.radialBounds = [1 15]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust for each subject [M P blank]
task.unitSigma = [.2 .4 .50];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

