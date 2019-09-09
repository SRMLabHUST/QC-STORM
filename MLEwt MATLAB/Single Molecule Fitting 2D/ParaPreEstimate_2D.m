function oinf= ParaPreEstimate_2D(moleculeSub, ScalingCoef)


RegionLength=size(moleculeSub,1);


imgdat=double(moleculeSub);
backg=mean([imgdat(1,1:end) imgdat(end,1:end) imgdat(2:end-1,1)' imgdat(2:end-1,end)']);

Amplitude=max(imgdat(:))-backg;

X0=RegionLength/2;
Y0=RegionLength/2;
SigmaX = 1.3;

sigmax = 1/(2*SigmaX*SigmaX);

Amplitude = Amplitude/ScalingCoef.CoefA;
backg = backg/ScalingCoef.CoefB;
sigmax = sigmax/ScalingCoef.CoefS;

oinf=[Amplitude; X0; Y0; sigmax; backg];


