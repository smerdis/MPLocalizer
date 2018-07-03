function [p t task] = runMPLocalizer(subjectID, run)
% [p t task] = runMPLocalzer(subjectID, run)
%
% Main script for LGN M/P localizer.
%
% Based on rd_mpLocalizerColor_v4.m
%
% Rachel Denison

KbName('UnifyKeyNames') % allows referring to keys as 5% etc

if nargin==0
    subjectID = 'test';
    run = 1;
end

triggerOnKey = 1; % 1 to start with key, 0 to start with TTL pulse only

% ------------------------------------------------------------------------
% Experiment setup
% ------------------------------------------------------------------------
% Load experiment parameters
[p.Gen task] = mpLocalizerColorParamsGen;
nConds = length(p.Gen.condNames);
for iCond = 1:nConds
    p.(p.Gen.condNames{iCond}) = mpLocalizerColorParamsStim(p.Gen.condNames{iCond});
end

if p.Gen.saveFile==0
    fprintf('\n\nNOT SAVING the file from this run.\n\n')
end

% Load device numbers
switch p.Gen.testingLocation
    case 'laptop'
        devNums = findKeyboardDevNums;
        load('gamma.mat','gammaTable')
        fprintf('\n\nLoading gamma table and devNums ...\n\n')
    otherwise
        fprintf('\n\nWill not load gamma table or devNums.\n\n')
end

% ------------------------------------------------------------------------
% Screen setup
% ------------------------------------------------------------------------
% PTB-3 correctly installed and functional? Abort otherwise.
AssertOpenGL;

% Select screen with maximum id for output window:
screenid = max(Screen('Screens'));

% Open a fullscreen, onscreen window with some color background.
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
[win winRect] = PsychImaging('OpenWindow', screenid, 0); %128
HideCursor;

% Load the gamma table
if exist('gammaTable','var')
    Screen('LoadNormalizedGammaTable', win, gammaTable);
else
    fprintf('\n\nDid not find gamma table.\n\n')
end

% Colors
black = BlackIndex(win);  % Retrieves the CLUT color code for black.
white = WhiteIndex(win);  % Retrieves the CLUT color code for white.
gray = (black + white) / 2;  % Computes the CLUT color code for gray.
bgColor = gray;

% Enable alpha-blending, set it to desired blend equation.
% Screen('BlendFunction', win, GL_ONE, GL_ONE); % if alpha blending the gabors
Screen('BlendFunction',win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % if using a transparent mask

% Query frame duration, define refresh rate
ifi = Screen('GetFlipInterval', win);
refrate = 1/ifi

% ------------------------------------------------------------------------
% Initial stimulus setup
% ------------------------------------------------------------------------
% Retrieve size of window in pixels
[w, h] = RectSize(winRect);
p.Gen.screenRes = [h w];
fprintf('\nSetting screen size to %i x %i.\n', h, w);

p.Gen.pixelsPerDegree = ang2pix(1.0, p.Gen.screenSize(1), p.Gen.screenRes(1), ...
    p.Gen.viewDist, 'central');

% Make polar grids used to specify stimulus layout
[r theta p.Gen] = makePolarMeshGrids(p.Gen);

% Make fixation
fix = ones(size(r))*bgColor;
fix(r < p.Gen.fixSize) = white;
[fixx fixy] = find(fix==white);
fixLBounds = min([fixx fixy]);
fixUBounds = max([fixx fixy]);
fixPatch = fix(fixLBounds(1):fixUBounds(1), fixLBounds(2):fixUBounds(2));
fixtex = Screen('MakeTexture', win, fixPatch);
blackfixtex = Screen('MakeTexture', win, white - fixPatch);

% ------------------------------------------------------------------------
% Define gabor properties and layout
% ------------------------------------------------------------------------
if isempty(p.Gen.gratingSize)
    gratingSize = pix2ang([h w], p.Gen.screenSize(1), ...
        p.Gen.screenRes(1), p.Gen.viewDist, 'central');
    p.Gen.gratingSize = gratingSize;
else
    gratingSize = p.Gen.gratingSize;
end

for iCond = 1:nConds
    cond = p.Gen.condNames{iCond};
    
    % Initial parameters of gabors:
    % Frequency of sine grating:
    freq(iCond) = p.(cond).spatialFrequency;
    % Color type:
    colorType{iCond} = p.(cond).colorType;
    % Contrast of grating:
    contrast(iCond) = p.(cond).contrast;
    
    redWeight(iCond) = p.(cond).redWeight;
    greenWeight(iCond) = p.(cond).greenWeight;
end

% ------------------------------------------------------------------------
% Set up stimulus timing
% ------------------------------------------------------------------------
for iCond = 1:nConds
    cond = p.Gen.condNames{iCond};
    
    % Wedge duration
    t.Gen.wedgeDuration = p.Gen.cycleDuration/p.Gen.nWedgePhases;

    % Wedge phase shift request times (within run (s))
    t.Gen.wedgeRequests = 0:(t.Gen.wedgeDuration+p.Gen.responseDuration):...
        (t.Gen.wedgeDuration+p.Gen.responseDuration)*(p.Gen.numCycles-1);
    
    fprintf('\nRequested flick rate for condition %s: %i Hz\n', ...
        cond, p.(cond).flickRate);
    
    % Number of flicks between orientation changes
    t.(cond).nFlicksPerOrientation = round(p.(cond).flickRate/p.(cond).orientationChangeRate);
    
    % Flick request times (within wedge (s))
    t.(cond).flickDuration = 1/p.(cond).flickRate;
    t.(cond).flickRequests = 0:t.(cond).flickDuration:(t.Gen.wedgeDuration-t.(cond).flickDuration);
    
    % Constrast step request times (within flick (s))
    t.(cond).cStepDuration = 1/p.(cond).cStepRate;
    t.(cond).cStepRequests = 0:t.(cond).cStepDuration:(t.(cond).flickDuration-t.(cond).cStepDuration);
    
    % Defining the contrast steps
    t.(cond).nCStepsPerFlick = t.(cond).flickDuration/t.(cond).cStepDuration; % must be even!
    if t.(cond).nCStepsPerFlick==1
        p.(cond).cSteps = 1;
    elseif mod(t.(cond).nCStepsPerFlick,2)
        error('t.nCStepsPerFlick must be even or 1!')
    else
        switch p.Gen.temporalWaveform
            case 'triangular'
                p.(cond).cStepSize = 1/(t.(cond).nCStepsPerFlick/2);
                cStepsAscending = 0:p.(cond).cStepSize:1;
                cStepsDescending = 1:-p.(cond).cStepSize:0;
                p.(cond).cSteps = [cStepsAscending(1:end-1) cStepsDescending(1:end-1)];
            case 'sin'
                grid = 0:pi/t.(cond).nCStepsPerFlick:pi;
                p.(cond).cSteps = sin(grid(1:end-1));
            otherwise
                error('p.Gen.temporalWaveform not found')
        end
    end
end

fprintf('\nRequested block times:\n')
disp(t.Gen.wedgeRequests + p.Gen.blankDuration(1));

% ------------------------------------------------------------------------
% Set up gabor textures and mask texture
% ------------------------------------------------------------------------
fprintf('\nGenerating textures ...\n')
phases = [0 pi];
orientations = p.Gen.gratingOrientations;
nOrientations = numel(orientations);

for iCond = 1:nConds
    cond = p.Gen.condNames{iCond};
    for iPhase = 1:numel(phases)
        phase = phases(iPhase);
        for iOrient = 1:nOrientations
            orientation = orientations(iOrient);
            for iFrame = 1:numel(p.(cond).cSteps)
                cStep = p.(cond).cSteps(iFrame);
                gratingIm = buildColorGrating(p.Gen.pixelsPerDegree, gratingSize, ...
                    freq(iCond), orientation, phase, ...
                    cStep*p.(cond).contrast, 0, colorType{iCond}, ...
                    redWeight(iCond), greenWeight(iCond));
                gratingtex{iCond}(iFrame,iPhase,iOrient) = Screen('MakeTexture', win, gratingIm*white);
            end
        end
    end
end

p.Gen.gratingPx = size(gratingIm);

% ------------------------------------------------------------------------
% Create mask
% ------------------------------------------------------------------------
fprintf('\nGenerating mask ...\n')
% Create wedge mask and/or masked centers image
switch p.Gen.maskType
    case 'transparentWedge'
        % Mask using transparency (2-layer mask)
        maskw = wedgify({zeros(h,w)}, r, theta, p.Gen, 0.5, 1);
        maskwAlphaLayer = ones(size(maskw{1}));
        maskwAlphaLayer(maskw{1}==0) = 0; % 0 is fully transparent
        
        mask(:,:,1) = maskw{1}*white;
        mask(:,:,2) = maskwAlphaLayer*white;
        
    case 'transparentGaussian'
        % Make transparent Gaussian mask, 2-layer alpha
        sigmaPx = p.Gen.gratingPx.*p.Gen.maskSigmaProp;
        maskGauss = make2DGaussianOval(w, h, w/2, h/2, ...
            sigmaPx(2), sigmaPx(1), p.Gen.maskAmp);
        maskGauss(maskGauss>1) = 1;
        
        mask(:,:,1) = ones(h, w)*bgColor;
        mask(:,:,2) = (1-maskGauss)*white;
        
    case 'rectangleBlur'
        recPx = p.Gen.gratingPx*p.Gen.apertureProp;
        rec = ones(h,w);
        rec(round((h-recPx(1))/2+1:(h-recPx(1))/2+recPx(1)),...
            round((w-recPx(2))/2+1:(w-recPx(2))/2+recPx(2))) = 0;
        mask(:,:,1) = ones(h, w)*bgColor;
        mask(:,:,2) = blurMaskEdges(rec,p.Gen.maskSmoothSize)*white;
        
    otherwise
        error('Mask type not found.')
end

% Make mask texture
masktex = Screen('MakeTexture', win, mask);

% ------------------------------------------------------------------------
% Set up target timing
% ------------------------------------------------------------------------
nTargetsOptions = 0:task.maxTargets;
for iBlock = 1:p.Gen.numCycles
    % randomly select the number of targets for this block
    nTargetsIdxs = randperm(task.maxTargets+1);
    task.nTargets(iBlock) = nTargetsOptions(nTargetsIdxs(1));
    targetPeriodDur = p.Gen.cycleDuration-task.targetCushion*2;
    
    % loop until there are no overlapping targets
    overlap = ones(1,task.maxTargets-1);
    while any(overlap)
        targetISIs = abs(normrnd(targetPeriodDur/(task.maxTargets+1),...
            task.ISISuggestedSigma,[task.maxTargets+1 1]));
        cumTargetISIs = cumsum(targetISIs);
        targetTimes = cumTargetISIs(1:task.maxTargets);
        targetTimesScaled = task.targetCushion + targetTimes*(targetPeriodDur/cumTargetISIs(end));
        for iTarget = 1:task.maxTargets-1
            overlap(iTarget) = targetTimesScaled(iTarget+1) < ...
                targetTimesScaled(iTarget) + task.targetCushion;
        end
    end
    task.targetTimesScaled{iBlock} = targetTimesScaled;
    
    % loop until the right number of targets are selected by binornd
    nTargetsNow = NaN;
    while nTargetsNow ~= task.nTargets(iBlock)
        task.targetOn{iBlock} = logical(binornd(1,.5,[task.maxTargets 1]));
        nTargetsNow = nnz(task.targetOn{iBlock});
    end
    
    task.targetPresentationTimes{iBlock} = ...
        task.targetTimesScaled{iBlock}(task.targetOn{iBlock});
    task.correctKey(iBlock) = ...
        p.Gen.keyCodes(task.nTargets(iBlock)==task.responseOptions);
    
    t.Gen.targetOnRequests{iBlock} = task.targetPresentationTimes{iBlock};
    t.Gen.targetOffRequests{iBlock} = t.Gen.targetOnRequests{iBlock} + task.targetDur;
end

fprintf('\nNumber of targets per block:\n')
disp(task.nTargets);

% ------------------------------------------------------------------------
% Set up target textures
% ------------------------------------------------------------------------
for iBlock = 1:p.Gen.numCycles
    
    condIdx = p.Gen.condOrder(iBlock);
    cond = p.Gen.condNames{condIdx};
    
    for iTarget = 1:task.nTargets(iBlock)
        [tr tth tx ty] = chooseTargetPosition(r, theta, task.radialBounds);
        task.targetPos{iBlock,iTarget} = [tx ty];
        unitSigma = task.unitSigma(p.Gen.condOrder(iBlock));
        task.targetSigma{iBlock,iTarget} = setTargetSigma(tr, unitSigma, p.Gen);
        
        targetImage(:,:,1) = ones(h, w)*p.(cond).targetColor(1);
        targetImage(:,:,2) = ones(h, w)*p.(cond).targetColor(2);
        targetImage(:,:,3) = ones(h, w)*p.(cond).targetColor(3);
        targetImage(:,:,4) = make2DGaussian(w, h, tx, ty, ...
            task.targetSigma{iBlock,iTarget}, task.targetAmp); % transparent layer, 0 is completely transparent
        
        targettex{iBlock,iTarget} = Screen('MakeTexture', win, targetImage*white);
    end
end

% ------------------------------------------------------------------------
% Preliminary save
% ------------------------------------------------------------------------
save('tempdata.mat', 'p', 't', 'task', 'subjectID', 'run')

% ------------------------------------------------------------------------
% Subject and scanner ready?
% ------------------------------------------------------------------------
% Fill background with gray
Screen('FillRect', win, bgColor);

% Press any key to begin
readyMessage = 'READY?\n\nPress any button to begin.';
DrawFormattedText(win, readyMessage, 'center', 'center');
Screen('Flip', win);
switch p.Gen.testingLocation
    case 'laptop'
        KbWait(devNums.Keypad);
    otherwise
        KbWait;
end

% Wait for scanner trigger (KbWait checks every 5 ms)
Screen('DrawTexture', win, fixtex);
DrawFormattedText(win, 'Waiting for scanner ...','center', p.Gen.screenRes(1)/2-50);
Screen('Flip', win);

if triggerOnKey
    % FOR KEY START
    switch p.Gen.testingLocation
        case 'laptop'
            t.Gen.trigger = KbWait(devNums.Keypad);
        otherwise
            t.Gen.trigger = KbWait;
    end
else
    % FOR SCANNER START
    keyPressed = '0';
    while ~strcmp(keyPressed,'5')
        [keyIsDown secs keyCode] = PsychHID('KbCheck',devNums.TTLPulse);
        keyPressed = KbName(keyCode);
        if length(keyPressed) > 1
            keyPressed = keyPressed(1);
        end
    end
    t.Gen.trigger = secs;
end

if p.Gen.blankDuration(1) > 0
    % Show fixation
    Screen('DrawTexture', win, fixtex);
    Screen('Flip', win);
end

% ------------------------------------------------------------------------
% Display stimulus
% ------------------------------------------------------------------------
count = 0;
wedgeIdx = 1;
t.Gen.flipsHeaders = {'vbl','iCond','wedgeIdx','flickIdx','cStepIdx','targetIsOn'};

while wedgeIdx <= length(t.Gen.wedgeRequests)

    % stimulus initializations
    iCond = p.Gen.condOrder(wedgeIdx);
    cond = p.Gen.condNames{iCond};
    wedgeIdx = wedgeIdx + 1;
    flickIdx = 1;
    rotIdx = 1;
    rotAngleIdxs = randperm(nOrientations);
    p.Gen.rotAngles{wedgeIdx-1} = p.Gen.gratingOrientations(rotAngleIdxs);
    imagetex = gratingtex{iCond}(:,mod(flickIdx,2)+1,rotAngleIdxs(rotIdx));
    
    % task initializations
    targetOnFlick = zeros(length(t.(cond).flickRequests), task.maxTargets);
    targetOnCStep = cell(1,length(t.(cond).flickRequests));
    
    for iTarget = 1:task.nTargets(wedgeIdx-1)        
        targetFirstFlick = find(diff(t.(cond).flickRequests < ...
            t.Gen.targetOnRequests{wedgeIdx-1}(iTarget)));
        targetLastFlick = find(diff(t.(cond).flickRequests < ...
            t.Gen.targetOffRequests{wedgeIdx-1}(iTarget)));
        targetOnFlick(targetFirstFlick:targetLastFlick,iTarget) = 1;
        targetOnFlickIdxs = find(targetOnFlick(:,iTarget));
        
        for iFlick = 1:length(targetOnFlickIdxs)         
            flick = targetOnFlickIdxs(iFlick);
            targetOnCStep{flick}(:,iTarget) = zeros(length(t.(cond).cStepRequests),1);
            
            targetFirstCStep = find(diff(t.(cond).flickRequests(flick) + ...
                t.(cond).cStepRequests < t.Gen.targetOnRequests{wedgeIdx-1}(iTarget)));
            if isempty(targetFirstCStep)
                targetFirstCStep = 1;
            end
            targetLastCStep = find(diff(t.(cond).flickRequests(flick) + ...
                t.(cond).cStepRequests < t.Gen.targetOffRequests{wedgeIdx-1}(iTarget)));
            if isempty(targetLastCStep)
                targetLastCStep = length(t.(cond).cStepRequests);
            end
            targetOnCStep{flick}(targetFirstCStep:targetLastCStep,iTarget) = 1;
        end
    end
    % store target flicks and cSteps
    task.targetOnFlick{wedgeIdx-1} = targetOnFlick;
    task.targetOnCStep{wedgeIdx-1} = targetOnCStep;

    while flickIdx <= length(t.(cond).flickRequests)

        if mod(flickIdx, t.(cond).nFlicksPerOrientation)==0
            % Compute new random orientation for each patch in next frame:
            rotIdx = rotIdx + 1;
        else
            % Increment phase-shift of each gabor by 180 deg. per flick
            imagetex = gratingtex{iCond}(:,mod(flickIdx,2)+1,rotAngleIdxs(rotIdx));
        end

        flickIdx = flickIdx + 1;
        cStepIdx = 1;

        while cStepIdx <= length(t.(cond).cStepRequests)

            cStepIdx = cStepIdx + 1;

            % Draw grating
            Screen('DrawTexture', win, imagetex(cStepIdx-1));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % task code
            % Draw target, if present (timing resolved at level of flicks and cSteps)            
            if any(targetOnFlick(flickIdx-1,:),2) && ...
                    any(targetOnCStep{flickIdx-1}(cStepIdx-1,:),2)
                Screen('DrawTexture', win, targettex{wedgeIdx-1, find(targetOnFlick(flickIdx-1,:))});
                targetIsOn = 1;
            else
                targetIsOn = 0;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Draw transparent mask
            Screen('DrawTexture', win, masktex);

            % Draw fixation
            Screen('DrawTexture', win, fixtex);

            % Mark drawing ops as finished, so the GPU can do its drawing job while
            % we can compute updated parameters for next animation frame.
            Screen('DrawingFinished', win);


            % Flip screen when time for next contrast step
            vbl = Screen('Flip', win, t.Gen.trigger + p.Gen.blankDuration(1) + ...
                t.Gen.wedgeRequests(wedgeIdx-1) + t.(cond).flickRequests(flickIdx-1) + ...
                t.(cond).cStepRequests(cStepIdx-1) - ifi/2);

            % Next loop iteration...
            count = count + 1;

            % Store time of screen flip, the flick idx, and the wedge idx on that flick
            t.Gen.flips(count,:) = [vbl iCond wedgeIdx-1 flickIdx-1 cStepIdx-1 targetIsOn];

        end
    end
    
    % Collect response 
    Screen('FillRect', win, bgColor);
    Screen('DrawTexture', win, blackfixtex);
    vbl = Screen('Flip', win, t.Gen.trigger + p.Gen.blankDuration(1) + ...
        t.Gen.wedgeRequests(wedgeIdx-1) + t.Gen.wedgeDuration - ifi/2);
    
    RT = NaN;
    keyPressed = NaN;
    while GetSecs < t.Gen.trigger + p.Gen.blankDuration(1) + ...
                t.Gen.wedgeRequests(wedgeIdx-1) + t.Gen.wedgeDuration ...
                + p.Gen.responseDuration - p.Gen.trialCushion - ifi/2
        %switch p.Gen.testingLocation
        %    case 'laptop'
        %        [keyIsDown secs keyCode] = PsychHID('KbCheck',devNums.Keypad);
        %    otherwise
        [keyIsDown secs keyCode] = KbCheck;
        %end
%         WaitSecs(.01);
        if keyIsDown 
            RT = secs - vbl;
            keyPressed = find(keyCode);
            if length(keyPressed) > 1
                keyPressed = keyPressed(1);
            end
        end
    end
    task.RTs(wedgeIdx-1) = RT;
    task.responseKey(wedgeIdx-1) = keyPressed;
    task.acc(wedgeIdx-1) = ...
        keyPressed==task.correctKey(wedgeIdx-1);
    
    Screen('DrawTexture', win, fixtex);
    vbl = Screen('Flip', win, t.Gen.trigger + p.Gen.blankDuration(1) + ...
        t.Gen.wedgeRequests(wedgeIdx-1) + t.Gen.wedgeDuration + ...
        p.Gen.responseDuration - p.Gen.trialCushion);
end

% ------------------------------------------------------------------------
% Display final blank period
% ------------------------------------------------------------------------
if p.Gen.blankDuration(2) > 0
    Screen('FillRect', win, bgColor);
    Screen('DrawTexture', win, fixtex);

    % Show fixation and take timestamp of stimulus offset
    vbl = Screen('Flip', win, t.Gen.trigger + p.Gen.blankDuration(1) + ...
        (t.Gen.wedgeDuration+p.Gen.responseDuration)*p.Gen.numCycles - ifi/2);
    t.Gen.endStim = vbl;

    % End of final blank period
    vbl = Screen('Flip', win, t.Gen.trigger + p.Gen.blankDuration(1) + ...
        (t.Gen.wedgeDuration+p.Gen.responseDuration)*p.Gen.numCycles + p.Gen.blankDuration(2) - ifi/2);
    t.Gen.endRun = vbl;
else
    % Last flip to take end timestamp and for stimulus offset
    vbl = Screen('Flip', win, t.Gen.trigger + p.Gen.blankDuration(1) + ...
        (t.Gen.wedgeDuration+p.Gen.responseDuration)*p.Gen.numCycles - ifi/2);
    t.Gen.endStim = vbl;
    t.Gen.endRun = vbl;
end

% ------------------------------------------------------------------------
% Timing stats
% ------------------------------------------------------------------------
count
avgfps = count / (t.Gen.endStim - t.Gen.trigger - p.Gen.blankDuration(1))

t.Gen.count = count;
t.Gen.avgfps = avgfps;

for iBlock = 1:wedgeIdx-1
    blockStartIdx = find(t.Gen.flips(:,3)==iBlock, 1, 'first');
    blockEndIdx = find(t.Gen.flips(:,3)==iBlock, 1, 'last');
    blockDuration = t.Gen.flips(blockEndIdx,1) - t.Gen.flips(blockStartIdx,1);
    nFlips = blockEndIdx - blockStartIdx;
    t.Gen.avgfpsByBlock(iBlock) = blockDuration/nFlips;
end

fprintf('\nAverage fps, by block: [condition fps]\n\n')
disp([p.Gen.condOrder' t.Gen.avgfpsByBlock'])

t.Gen.flips = [t.Gen.flips(:,1)-t.Gen.trigger t.Gen.flips(:,2:end)]; % relative to trigger
t.Gen.flipsRel = [[0 99 99 99 99 99]; t.Gen.flips; [t.Gen.endStim-t.Gen.trigger 99 99 99 99 99]]; % first event represents the trigger
t.Gen.flicksRel = t.Gen.flipsRel(find(diff(t.Gen.flipsRel(:,4)))+1, [1 4]);
t.Gen.wedgesRel = t.Gen.flipsRel(find(diff(t.Gen.flipsRel(:,3)))+1, (1:3));

fprintf('\nActual block times + end time:\n\n')
disp(t.Gen.wedgesRel)
disp(t.Gen.endRun - t.Gen.trigger)

figure
plot(diff(t.Gen.flicksRel(:,1)),'.')
title('Flick intervals')

% ------------------------------------------------------------------------
% Behavioral performance
% ------------------------------------------------------------------------
fprintf('\nBehavioral performance:\n\n')
for iCond = 1:length(p.Gen.condNames)
    trialAcc = task.acc(p.Gen.condOrder==iCond);
    task.nTrialsInCond(iCond) = numel(trialAcc);
    task.condAcc(iCond) = mean(trialAcc);
    fprintf('%s \t%.2f   n=%d\n', ...
        p.Gen.condNames{iCond}, task.condAcc(iCond), ...
        task.nTrialsInCond(iCond))
end 
fprintf('\n')

% ------------------------------------------------------------------------
% Behavioral feedback
% ------------------------------------------------------------------------
% Overall accuracy for all conditions but blank
overallAcc = ...
    mean(task.acc(p.Gen.condOrder~=find(strcmp(p.Gen.condNames,'blank'))));
feedbackText = sprintf('Accuracy: %d%%', round(overallAcc*100));
DrawFormattedText(win, feedbackText, 'center', 'center')
Screen('Flip', win)
WaitSecs(3);

% ------------------------------------------------------------------------
% Close onscreen window
% ------------------------------------------------------------------------
Screen('CloseAll');
ShowCursor;

% ------------------------------------------------------------------------
% Clean up
% ------------------------------------------------------------------------
% Timestamp, subject ID
p.Gen.whenSaved = datestr(now);
p.Gen.subjectID = subjectID;
p.Gen.run = run;

% Save file
if p.Gen.saveFile
    saveFile = sprintf('data/mpLocalizer_%s_run%02d_%s.mat', subjectID, run, datestr(now,'yyyymmdd'));
    if exist(saveFile, 'file')
        fprintf('\n\nFile already exists. Saving generic-named file ...\n\n')
        save(sprintf('DATA_%s.mat',datestr(now,'yyyymmdd')), 'p', 't', 'task', 'subjectID', 'run')
    else
        if ~exist('./data','dir')
            mkdir('data')
        end
        save(saveFile, 'p', 't', 'task', 'subjectID', 'run')
        fprintf('\n\nFile saved.\n\n')
    end

    % write 3-col events files
    events_fn = sprintf('data/sub-%s_ses-%s_task-mp_run-%02d_events.tsv', subjectID, datestr(now,'yyyymmdd'), run);
    blocks_in_order = p.Gen.condNames(p.Gen.condOrder);
    events_file_contents = 'onset\tduration\ttrial_type\n' ;
    if p.Gen.cycleDuration == 18 %3T
        time_between_onsets = 20.25 ;
    elseif p.Gen.cycleDuration == 16 % 7T
        time_between_onsets = 18 ;
    else
        fprintf('\nERROR: cycleDuration is not 16 or 18! Assuming 18...\n') ;
        time_between_onsets = 20.25 ;
    end
    for i=1:length(blocks_in_order)
        events_file_contents = [events_file_contents sprintf('%.02f\t%d\t%s\n',(i-1)*time_between_onsets, p.Gen.cycleDuration, blocks_in_order{i})] ;
    end
    events_fid = fopen(events_fn, 'w');
    fprintf(events_fid, events_file_contents);
    fclose(events_fid);
end

% Done.
return;

