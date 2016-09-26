function [p t] = rd_hemifieldMapping_v3(subjectID, run, steadyState)
%
% function [p t] = rd_hemifieldMapping(subjectID, run, [steadyState])
% 
% Left/right hemifield localizer, block design.
%
% INPUTS:
% subjectID is a string with a subject identifier
% run is the run number
% steadyState (OPTIONAL) should be 'steadyState' if you want to do a long,
% full-field block.
%
% Modified from rd_hemifieldCheckerboardEvents_v2.m
%
% v3 implements custom height mask
%
% Rachel Denison
% 22 July 2012

% triggerOnKey = 1; % use keyboard press as the trigger
triggerOnKey = [];
while isempty(triggerOnKey) || ~ismember(triggerOnKey, [1 0])
    triggerOnKey = input('\nUse key press as trigger? (1 for yes, 0 for TTL only) ');
end

% ------------------------------------------------------------------------
% Experiment setup
% ------------------------------------------------------------------------
% Load experiment parameters
p = hemifieldParams;

if exist('steadyState','var') && strcmp(steadyState,'steadyState')
    p = hemiUpdateSteadyStateParams(p);
end 

% Check params
fprintf('\nTesting location: %s\n', p.testingLocation)
fprintf('Screen size: [%.02f %.02f]\n', p.screenSize)
fprintf('Viewing distance: %.02f\n\n', p.viewDist)
ok = input('Ok to go on? (n if not): ', 's');
if strcmp(ok,'n')
    error('Exiting ... check hemifieldParams.m')
end

if p.saveFile==0
    fprintf('\n\nNOT SAVING the file from this run.\n\n')
end

% Load device numbers
switch p.testingLocation
    case '7T'
        load('7TASGamma15b.mat','gammaTableMeanBottom')
        gammaTable = gammaTableMeanBottom;
        devNums = findKeyboardDevNums_7T;
        if isempty(devNums.ButtonBox)
            error('Could not find button box! Please check findKeyboardDevNums.')
        end
        if isempty(devNums.TTLPulse)
            error('Could not find TTL device! Please check findKeyboardDevNums.')
        end
        fprintf('\n\nLoading 7T-AS gamma table and devNums ...\n\n')
    case '3T'
        load('BIC3TAvotecGamma.mat','gammaTable')
        devNums = findKeyboardDevNums_scanner;
        if isempty(devNums.ButtonBox)
            error('Could not find button box! Please check findKeyboardDevNums.')
        end
        if isempty(devNums.TTLPulse)
            error('Could not find TTL device! Please check findKeyboardDevNums.')
        end
    otherwise
        fprintf('\n\nNot testing in scanner ... will not load gamma table or devNums.\n\n')
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
% Stimulus setup
% ------------------------------------------------------------------------
% Retrieve size of window in pixels
[w, h] = RectSize(winRect);
p.screenRes = [h w];
fprintf('\nSetting screen size to %i x %i.\n', h, w);

% Make polar grids used to specify stimulus layout
[r theta p] = makePolarMeshGrids(p);

% Make radial checkerboards (contrast reversed)
[img{1} img{2}] = makeRadialCheckerboards(p, r, theta);
for iFlick = 1:length(img)
    imagetex{iFlick} = Screen('MakeTexture', win, img{iFlick}*white);
end

% Make fixation
fix = ones(size(r))*bgColor;
fix(r < p.fixSize) = white;
[fixx fixy] = find(fix==white);
fixLBounds = min([fixx fixy]);
fixUBounds = max([fixx fixy]);
fixPatch = fix(fixLBounds(1):fixUBounds(1), fixLBounds(2):fixUBounds(2));
fixtex = Screen('MakeTexture', win, fixPatch);

% Make mask images
maskw = wedgify({zeros(size(img{1}))}, r, theta, p, 0.5, 1);

% Make mask using transparency (2-layer mask)
for iMask = 1:length(maskw)
    maskwa{iMask}(:,:,1) = maskw{iMask};
    maskwAlphaLayer = ones(size(maskw{iMask}));
    maskwAlphaLayer(maskw{iMask}==0) = 0; % 0 is fully transparent
    maskwa{iMask}(:,:,2) = maskwAlphaLayer;
    masktex{iMask} = Screen('MakeTexture', win, maskwa{iMask}*white);
end

% Make custom height mask if requested (2-layer mask)
if ~isempty(p.customMaskHeight)
    maskHeightPx = ang2pix(p.customMaskHeight,...
        p.screenSize(2), p.screenRes(2), p.viewDist, 'central');
    p.customMask(:,:,1) = ones(h,w)*bgColor;
    p.customMask(:,:,2) = ones(h,w)*white;
    p.customMask(round(h/2-maskHeightPx/2):round(h/2+maskHeightPx/2),:,2) = 0; % transparent
    customMaskTex = Screen('MakeTexture', win, p.customMask);
end


% ------------------------------------------------------------------------
% Set up timing
% ------------------------------------------------------------------------
% Wedge duration
t.wedgeDuration = p.cycleDuration/p.nWedgePhases;

% Wedge phase shift request times (within run (s))
t.wedgeRequests = 0:t.wedgeDuration:t.wedgeDuration*...
    (p.numCycles*p.nWedgePhases-1);

% Flick request times (within wedge (s))
t.flickDuration = 1/p.flickRate;
t.flickRequests = 0:t.flickDuration:(t.wedgeDuration-t.flickDuration);

% ------------------------------------------------------------------------
% Preliminary save
% ------------------------------------------------------------------------
save('tempdata.mat', 'p', 't', 'subjectID', 'run')

% ------------------------------------------------------------------------
% Subject and scanner ready?
% ------------------------------------------------------------------------
% Fill background with gray
Screen('FillRect', win, bgColor);

% Press any key to begin
readyMessage = 'READY?\n\nPress any button to begin.';
DrawFormattedText(win, readyMessage, 'center', 'center');
Screen('Flip', win);
switch p.testingLocation
    case {'3T','7T'}
        KbWait(devNums.ButtonBox);
    otherwise
        KbWait;
end

% Wait for scanner trigger (KbWait checks every 5 ms)
Screen('DrawTexture', win, fixtex);
DrawFormattedText(win, 'Waiting for scanner ...','center', p.screenRes(1)/2-50);
Screen('Flip', win);

if triggerOnKey
    % FOR KEY START
    [t.trigger, keyCode] = KbWait;
else
    % FOR SCANNER START
    keyPressed = [];
    while isempty(keyPressed) || keyPressed~=p.triggerCode
        [keyIsDown secs keyCode] = PsychHID('KbCheck',devNums.TTLPulse);
        keyPressed = find(keyCode);
        if length(keyPressed) > 1;
            keyPressed = keyPressed(1);
        end
    end
    t.trigger = secs;
end

% Draw initial blank screen
Screen('DrawTexture', win, fixtex);
t.startFlip = Screen('Flip', win);

% ------------------------------------------------------------------------
% Display stimulus
% ------------------------------------------------------------------------
count = 0;
wedgeIdx = 0;

while wedgeIdx < length(t.wedgeRequests)
    wedgeIdx = wedgeIdx + 1;
    flickIdx = 0;

    while flickIdx < length(t.flickRequests)
        flickIdx = flickIdx + 1;

        % Draw image
        Screen('DrawTexture', win, imagetex{2-mod(flickIdx,2)})
        
        % Draw transparent mask
        Screen('DrawTexture', win, masktex{2-mod(wedgeIdx,2)});
        
        % Draw custom mask if requested
        if ~isempty(p.customMaskHeight)
            Screen('DrawTexture', win, customMaskTex);
        end
        
        % Flip screen
        vbl = Screen('Flip', win, t.trigger + ...
            t.wedgeRequests(wedgeIdx) + ...
            t.flickRequests(flickIdx) - ifi/2);
        
        % Next loop iteration...
        count = count + 1;
        
        % Store time of screen flip, the flick idx, and the wedge idx on that flick
        t.flips(count,:) = [vbl wedgeIdx flickIdx];
    end
end

% Final flip
Screen('FillRect', win, bgColor);
Screen('DrawTexture', win, fixtex);
t.endFlip = Screen('Flip', win, t.trigger + ...
    p.cycleDuration*p.numCycles - ifi/2);
t.flipsHeaders = {'vbl','wedgeIdx','flickIdx'};

% ------------------------------------------------------------------------
% Close onscreen window
% ------------------------------------------------------------------------
Screen('CloseAll');
ShowCursor;

% ------------------------------------------------------------------------
% Stimulus timing
% ------------------------------------------------------------------------
t.flipsRel = t.flips(:,1) - t.trigger;
t.flipsRel(:,2:3) = t.flips(:,2:3);
t.flicksRel = t.flipsRel(find(diff(t.flipsRel(:,3)))+1, [1 3]);
t.wedgesRel = t.flipsRel(find(diff(t.flipsRel(:,2)))+1, (1:2));

t.flickDurs = diff([0 t.flicksRel(:,1)' t.endFlip-t.trigger]);
t.wedgeDurs = diff([0 t.wedgesRel(:,1)' t.endFlip-t.trigger]);

fprintf('\nActual wedge durations:\n\n')
disp(t.wedgeDurs)

figure
plot(diff(t.flicksRel(:,1)),'.')
title('Flick intervals')

% ------------------------------------------------------------------------
% Clean up
% ------------------------------------------------------------------------
% Timestamp, subject ID
p.whenSaved = datestr(now);
p.subjectID = subjectID;
p.run = run;

% Save file
if p.saveFile
    saveFile = sprintf('data/%s_run%02d_HemifieldMapping_%s.mat', subjectID, run, datestr(now,'yyyymmdd'));
    if exist(saveFile, 'file')
        fprintf('\n\nFile already exists. Saving generic-named file ...\n\n')
        save(sprintf('DATA_%s.mat',datestr(now,'yyyymmdd')), 'p', 't', 'subjectID', 'run')
    else
        save(saveFile, 'p', 't', 'subjectID', 'run')
        fprintf('\n\nFile saved.\n\n')
    end
end

% Done.
return;


