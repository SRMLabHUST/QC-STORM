%% simu parameters
% ver1 use Gaussian PSF


TotalPhoton = 3000;
BackgroundPhoton = 450;

PSFSigmaX = 1.4017; % pixel
PSFSigmaY = 1.2815; % pixel
PixelSize = 137; % nm
    
QE = 0.95;
ReadNoise = 0;
TakeReadNoiseToBackground = 0;
IsEMCCD = 1;

ROISize = 15;

TotalPhoton=LocArry(:,7);
BackgroundPhoton=LocArry(:,8);
PSFSigmaX=LocArry(:,5);
PSFSigmaY=LocArry(:,6);

%%

ParaStd_all = GetLocPrec_CRLB3D(TotalPhoton, BackgroundPhoton, PSFSigmaX, PSFSigmaY, QE, ReadNoise, TakeReadNoiseToBackground, IsEMCCD, ROISize);

StdX0 = mean(ParaStd_all(:,2)) * PixelSize
StdY0 = mean(ParaStd_all(:,3)) * PixelSize
mean(ParaStd_all(:,4))* PixelSize
mean(ParaStd_all(:,5))* PixelSize
