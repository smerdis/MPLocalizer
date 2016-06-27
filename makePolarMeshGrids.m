function [r theta p] = makePolarMeshGrids(p)

[Mx My] = meshgrid(1:p.screenRes(2), 1:p.screenRes(1));

sz = size(Mx);

r = pix2ang(radialize(Mx,My), p.screenSize(1), p.screenRes(1), p.viewDist, 'radial'); % radial meshgrid

x0 = My - sz(1)/2;
y0 = Mx - sz(2)/2;
theta = atan2(x0,y0);
theta(theta<=0) = theta(theta<=0) + 2*pi;