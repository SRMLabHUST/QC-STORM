function fluo = Gaussian1( x0, y0, Sigma, ROISize)

HalfLen=floor(ROISize/2);

x=-HalfLen:HalfLen;
y=-HalfLen:HalfLen;



[X,Y]=meshgrid(x,y);

fluo = 1/(2*pi*Sigma^2)*exp(-((X-x0).^2 + (Y-y0).^2)/(2*Sigma^2));

