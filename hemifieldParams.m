function p = hemifieldParams

p.saveFile = 1;
p.testingLocation = '7T'; 

switch p.testingLocation
    case '7T'
        p.screenSize = [24 31]; % (cm) %% Measure!
        p.screenRes = [1200 1600];
        p.viewDist = 30; % (in) %% Measure!
        p.triggerCode = KbName('5%'); % keyCode for 5% is 34 
    case '3T' 
        p.screenSize = [7.17 9.25]; % (in) RD's measurement
        p.screenRes = [768 1024]; % 60 Hz. Must set screen to this!
        p.viewDist = 11.25; % (in) RD's measurement    
        % so screen is ? x ? deg visual angle
        p.triggerCode = KbName('5'); % keyCode for 5 is 93     
    case 'desk'
        p.screenSize = [12 14.75]; % (in)
        p.screenRes = [1024 1280];
        p.viewDist = 30; % (in)
    case 'laptop'
        p.screenSize = [9 13]; % (in)
        p.screenRes = [900 1440];
        p.viewDist = 21; % (in)
        p.triggerCode = KbName('5%'); % keyCode for 5% is 34
    otherwise
        error('Testing location not found!')
end

p.numCycles = 8.5; % 9 % a cycle means one left-right cycle
p.cycleDuration = 32; % (s) each left or right block is cycleDuration/2
p.nWedgePhases = 2; % hemifields
p.flickRate = 8.0;

p.wedgeRadial = [0.3 30]; % radius (internal and external) of wedge (deg) eg. [0.5 12]
p.wedgePolar = 180*(pi/180); % (rad) 
p.fixSize = 0.1; %0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% custom mask
p.customMaskHeight = 7; % in deg vis angle, [] if no custom mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.thetaCyc = 6; % Full cycles per pi (180 deg)
p.E = 0.05; 
p.A = 0.8; % E=0.5, A=0.3 based on figure 9 in Dumoulin and Wandell, 2007 or 0.05 and 0.8, based on Michael's measurements
p.contrast = 1;

p.mappingType = 'hemifieldLR';
p.maskType = 'transparentWedge';


