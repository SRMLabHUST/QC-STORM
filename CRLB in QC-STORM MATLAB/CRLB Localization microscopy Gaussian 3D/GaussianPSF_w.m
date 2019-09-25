function fluo = GaussianPSF_w(InputPara, ROISize)

TotalPn = InputPara(1);
x0 = InputPara(2);
y0 = InputPara(3);
SigmaX = InputPara(4);
SigmaY = InputPara(5);
Bkg   = InputPara(6);


fluo1 = Gaussian1(x0, y0, SigmaX, SigmaY, ROISize);
fluo1 = TotalPn*fluo1/sum(fluo1(:));

fluo = fluo1 + Bkg;

