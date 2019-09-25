%% simu parameters
% ver1 use Gaussian PSF

%%
PixelSize = 122; % nm

IsEMCCD = 1;
QE = 0.95;
ReadNoise = 1.0;
TakeReadNoiseToBackground = 1;

ROISize = 15;

%%

TotalFluoNum = size(LocArry, 1);

TotalPhoton = LocArry(:,7);
BackgroundPhoton = LocArry(:,8);

PSFSigmaX = LocArry(:,5); % pixel
PSFSigmaY = LocArry(:,6); % pixel


%%

CRLB_XY_nm = GetLocPrec_CRLB3D(TotalPhoton, BackgroundPhoton, PSFSigmaX, PSFSigmaY, PixelSize, QE, ReadNoise, TakeReadNoiseToBackground, IsEMCCD, ROISize);

pos = (CRLB_XY_nm(:,1)<40) & (CRLB_XY_nm(:,2)<40);

CRLB_XY_nm_1 = CRLB_XY_nm(pos,:);

save CRLB_XY_raw CRLB_XY_nm CRLB_XY_nm_1 IsEMCCD QE ReadNoise TakeReadNoiseToBackground ROISize PixelSize 

load('CalibrationData3D.mat');

Calib3d_fitresult = fitresult;
ZCrlb = GetZLocPrec_all(PSFSigmaX, PSFSigmaY, CRLB_XY_nm(:,1), CRLB_XY_nm(:,2), PixelSize, Calib3d_fitresult);
LocPrec_Z = mean(ZCrlb);

LocPrec_X = mean(CRLB_XY_nm_1(:,1));
LocPrec_Y = mean(CRLB_XY_nm_1(:,2));
LocPrec_Z = LocPrec_X + LocPrec_Y;

% save MeanLocPrec LocPrec_X LocPrec_Y LocPrec_Z
