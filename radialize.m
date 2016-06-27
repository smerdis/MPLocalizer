function r = radialize(Mx, My)

% function r = radialize(Mx, My)
%
% Takes X and Y meshgrid() outputs and returns a single radial grid.

sz = size(Mx);

r = sqrt( (My-sz(1)/2).^2 + (Mx-sz(2)/2).^2 );

return