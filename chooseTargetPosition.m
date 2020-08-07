function [r0 th0 x0 y0] = chooseTargetPosition(r, theta, radialBounds)
%
% function [r0 th0 x0 y0] = chooseTargetPosition(r, theta, radialBounds)
%
% r is a radial polar mesh grid in degrees vis ang
%
% theta is an angular polar mesh grid
%
% radialBounds is a two-element vector with inner and outer radial
%   boundaries in degrees visual angle

c = 0;

while c==0
    th0 = rand*2*pi; % polar angle
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
