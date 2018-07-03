function p = hemifieldParams

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

% addition by AM: include the TR and make sure block lengths align with it.
p.TR = 2.25; % for 3T using previous RD ep2d_neuro_... sequence

p.numCycles = 10.5; % a cycle means one left-right cycle, subtract .5 because reasons...
p.cycleDuration = 27; % (s) each left or right block is cycleDuration/2
assert(mod(p.numCycles*p.cycleDuration,p.TR)==0);
p.nWedgePhases = 2; % hemifields
p.flickRate = 8.0; % I (AM) believe this is half a flicker cycle, so 8.0 = 4Hz

p.wedgeRadial = [0.3 30]; % radius (internal and external) of wedge (deg) eg. [0.5 12]
p.wedgePolar = 180*(pi/180); % (rad) 
p.fixSize = 0.1; %0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% custom mask
p.customMaskHeight = []; % in deg vis angle, [] if no custom mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.thetaCyc = 6; % Full cycles per pi (180 deg)
p.E = 0.05; 
p.A = 0.8; % E=0.5, A=0.3 based on figure 9 in Dumoulin and Wandell, 2007 or 0.05 and 0.8, based on Michael's measurements
p.contrast = 1;

p.mappingType = 'hemifieldLR';
p.maskType = 'transparentWedge';


