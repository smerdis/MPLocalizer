function p = mpConnectivityColorParamsStim(stimType, Gen)

% Set r/g weights based on flicker photometry:
p.redWeight = 1;
p.greenWeight = 100/255;

% Note: for 60 Hz refrate, possible flicker freqs are [30 15 10 7.5 5 3 1.5 1]
p.stimType = stimType;

% Stim type M or P?
if any(regexp(p.stimType,'M')) % will catch leftM, rightM etc
    % M-type
    p.flickerFrequency = 15; % Hz
    p.spatialFrequency = 0.5; % cycles per degree
    p.colorType = 'bw';
    p.contrast = 1;
    p.cStepRate = 30;
    p.targetColor = [.5 .5 .5];
elseif any(regexp(p.stimType,'P')) % leftP, rightP
    % P-type
    p.flickerFrequency = 5; % Hz
    p.spatialFrequency = 2; % cycles per degree 
    p.colorType = 'rg';
    p.contrast = 1;
    p.cStepRate = 20;
    p.targetColor = [p.redWeight/2 p.greenWeight/2 0];
elseif any(regexp(p.stimType,'blank'))
    p.flickerFrequency = 0.5;
    p.spatialFrequency = 0.5;
    p.colorType = 'bw';
    p.contrast = 0;
    p.cStepRate = 1;
    p.targetColor = [.5 .5 .5];
else
    error('Stim type not found!')
end

p.flickRate = p.flickerFrequency*2; % flicks per s; 2 flicks make a full cycle
p.orientationChangeRate = length(Gen.gratingOrientations)/(Gen.cycleDuration) % orientation changes per s



