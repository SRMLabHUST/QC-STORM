function oinf= ParaPreEstimate_s3D(moleculeSub, ScalingCoef)

RegionLength = size(moleculeSub,1);


imgdat = double(moleculeSub);
Bkg = mean([imgdat(1,1:end) imgdat(end,1:end) imgdat(2:end-1,1)' imgdat(2:end-1,end)']);

Amplitude = max(imgdat(:)) - Bkg;

X0 = RegionLength/2;
Y0 = RegionLength/2;
SigmaX = 1.4;
SigmaY = 1.4;

SigmaX1 = 1/(2*SigmaX*SigmaX);
SigmaY1 = 1/(2*SigmaY*SigmaY);

Amplitude = Amplitude/ScalingCoef.CoefA;
Bkg = Bkg/ScalingCoef.CoefB;
SigmaX1 = SigmaX1/ScalingCoef.CoefS;
SigmaY1 = SigmaY1/ScalingCoef.CoefS;

oinf=[Amplitude; X0; Y0; SigmaX1; SigmaY1; Bkg];


