function fluo = Gaussian1(x0, y0, SigmaX, SigmaY, ROISize)

HalfLen = floor(ROISize/2);

x = -HalfLen:HalfLen;
y = -HalfLen:HalfLen;



[X, Y] = meshgrid(x, y);

fluo = 1/(2*pi*SigmaX*SigmaY)*exp(-((X - x0).^2/(2*SigmaX^2) + (Y - y0).^2/(2*SigmaY^2)));

