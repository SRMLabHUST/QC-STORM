function ModelSignal = EstimatedSignal_3D_2E(InInf, ScalingCoef, ROISize)

x=(1:ROISize)-0.5;
y=(1:ROISize)-0.5;

[X,Y] = meshgrid(x, y);


%%
A1 = InInf(1)*ScalingCoef.CoefA;
x1 = InInf(2);
y1 = InInf(3);

SigmaX1 = InInf(4)*ScalingCoef.CoefS;
SigmaY1 = InInf(5)*ScalingCoef.CoefS;

A2 = InInf(6)*ScalingCoef.CoefA;
x2 = InInf(7);
y2 = InInf(8);

SigmaX2 = InInf(9)*ScalingCoef.CoefS;
SigmaY2 = InInf(10)*ScalingCoef.CoefS;

Bkg = InInf(11)*ScalingCoef.CoefB;


Im1 = A1*exp(-((X-x1).*(X-x1)*SigmaX1 + (Y-y1).*(Y-y1))*SigmaY1);
Im2 = A2*exp(-((X-x2).*(X-x2)*SigmaX2 + (Y-y2).*(Y-y2))*SigmaY2);


ModelSignal = Im1 + Im2 + Bkg;
