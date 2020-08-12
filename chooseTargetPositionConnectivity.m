function [r0 th0 x0 y0] = chooseTargetPositionConnectivity(r, theta, radialBounds, cond)
%
% function [r0 th0 x0 y0] = chooseTargetPositionConnectivity(r, theta, radialBounds, cond)
%
% r is a radial polar mesh grid in degrees vis ang
%
% theta is an angular polar mesh grid
%
% radialBounds is a two-element vector with inner and outer radial
%   boundaries in degrees visual angle
%
% cond is the condition name, and will contain 'left' or 'right' to
% indicate which visual hemifield target should appear in.

c = 0;

while c==0
    if any(regexp(cond, 'left'))
        th0 = (rand*pi)+(pi/2); % polar angle
    elseif any(regexp(cond, 'right'))
        th0 = (rand*pi)-(pi/2); % polar angle
        if th0 < 0
            th0 = 2*pi+th0 ;
        end
    else
        th0 = rand*2*pi; % polar angle
    end
    r0 = rand*(radialBounds(2)-radialBounds(1)) + radialBounds(1); % ecc
    
    rmatch = find(round(r*10)==round(r0*10));
    thmatch = find(round(theta*30)==round(th0*30));
    
    matchim = ones(size(r));
    matchim(thmatch) = 0;
    matchim(rmatch) = 2;
    
    for i=1:length(rmatch)
        if any(rmatch(i)==thmatch)
            c = c+1;
            matchcoords(c) = rmatch(i);
        end
    end
%     fprintf('%d matches found\n', c)
end

matchim(matchcoords) = 10;
% imagesc(matchim)

[matchy matchx] = ind2sub(size(r),matchcoords);

% make sure matchx and matchy are even
if mod(numel(matchx),2)
    matchx = [matchx NaN];
end
if mod(numel(matchy),2)
    matchy = [matchy NaN];
end

x0 = matchx(round(numel(matchx))/2);
y0 = matchy(round(numel(matchy))/2);
