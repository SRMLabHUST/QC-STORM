function [fluo, SNR]  = GaussianPSF_w(InputPara, ROISize)

TotalPn   = InputPara(1);
x0    = InputPara(2);
y0    = InputPara(3);
Sigma = InputPara(4);
Bkg   = InputPara(5);

CenterROI=round(ROISize/2);


fluo1 = TotalPn*Gaussian1( x0, y0, Sigma, ROISize);

fluo = fluo1 + Bkg;

SNR=fluo1(CenterROI,CenterROI)/sqrt(fluo1(CenterROI,CenterROI)+Bkg);
