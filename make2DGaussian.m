function gaussian = make2DGaussian(w, h, x0, y0, sigma, gaussAmp)
%
% function gaussian = make2DGaussian(w, h, x0, y0, sigma, gaussAmp)


[x y] = meshgrid(1:w, 1:h);

gaussian = gaussAmp*exp(-((x-x0).^2/(2*sigma.^2))-((y-y0).^2/(2*sigma.^2)));