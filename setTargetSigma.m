function sigma = setTargetSigma(r0, unitSigma, p)
%
% function sigma = setTargetSigma(r0, unitSigma, p)
%
% r0 is the eccentricity of the target
% unitSigma is the standard deviation of the target in degrees visual angle
% p is a parameter structure with screen information

% convert sigma to pixelssca
unitSigmaPx = ...
    ang2pix(unitSigma, p.screenSize(1), p.screenRes(1), p.viewDist, 'central');

% linear scaling of target size with target eccentricity
sigma = unitSigmaPx*r0;