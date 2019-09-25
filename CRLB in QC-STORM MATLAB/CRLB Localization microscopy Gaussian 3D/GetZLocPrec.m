function ZStdError = GetZLocPrec(PSFSigmaX, PSFSigmaY, LocPrec_X, LocPrec_Y, PixelSize, Calib3d_fitresult)

d_SigmaX = LocPrec_X/PixelSize;
d_SigmaY = LocPrec_Y/PixelSize;

PSFSigmaX_simu = normrnd(PSFSigmaX,d_SigmaX,1,4000);
PSFSigmaY_simu = normrnd(PSFSigmaY,d_SigmaY,1,4000);

SigmaDiff = PSFSigmaX_simu.^2 - PSFSigmaY_simu.^2;

ZDepth_simu = 0.75*(Calib3d_fitresult.p1*SigmaDiff.^4 + Calib3d_fitresult.p2*SigmaDiff.^3 + Calib3d_fitresult.p3*SigmaDiff.^2 + Calib3d_fitresult.p4*SigmaDiff + Calib3d_fitresult.p5);

ZStdError = std(ZDepth_simu);