%% simu parameters
% ver1 use Gaussian PSF


TotalPhoton = 2000;
BackgroundPhoton = 70;

PSFSigmaX = 1.3; % pixel
PSFSigmaY = 1.3; % pixel
PixelSize = 100; % nm

QE = 0.72;
ReadNoise = 1.3;
TakeReadNoiseToBackground = 1;
IsEMCCD = 0;


%%

CRLB_X_nm = GetLocPrec_CRLB(TotalPhoton, BackgroundPhoton, PSFSigmaX, PixelSize, QE, ReadNoise, TakeReadNoiseToBackground, IsEMCCD)

