load('CameraInf.mat') % PixelSize 


LocPrec_X_pixel = 7/PixelSize;
LocPrec_Y_pixel = 7/PixelSize;

SigmaX = mean(LocArry(:,5));
SigmaY = mean(LocArry(:,6));

SigmaXVary = SigmaX + normrnd(0,LocPrec_X_pixel,1,10000);
SigmaYVary = SigmaY + normrnd(0,LocPrec_Y_pixel,1,10000);

SigmaDiff = SigmaXVary.^2 - SigmaYVary.^2;


ZDepthVary = 0.75*(fitresult.p1*SigmaDiff.^4 + fitresult.p2*SigmaDiff.^3 + fitresult.p3*SigmaDiff.^2 + fitresult.p4*SigmaDiff + fitresult.p5);

LocPrec_X_nm = std(SigmaXVary)*PixelSize
LocPrec_Y_nm = std(SigmaYVary)*PixelSize
LocPrec_Z_nm = std(ZDepthVary)

LocPrec_Z_nm/(LocPrec_X_nm+LocPrec_Y_nm)
