function ModelSignal = EstimatedSignal_2D_2E(InInf, ScalingCoef, ROISize)

x=(1:ROISize)-0.5;
y=(1:ROISize)-0.5;

[X,Y] = meshgrid(x, y);


%%
A1 = InInf(1)*ScalingCoef.CoefA;
x1 = InInf(2);
y1 = InInf(3);

A2 = InInf(4)*ScalingCoef.CoefA;
x2 = InInf(5);
y2 = InInf(6);

Sigma1 = InInf(7)*ScalingCoef.CoefS;
Bkg = InInf(8)*ScalingCoef.CoefB;


Im1 = A1*exp(-((X-x1).*(X-x1) + (Y-y1).*(Y-y1))*Sigma1);
Im2 = A2*exp(-((X-x2).*(X-x2) + (Y-y2).*(Y-y2))*Sigma1);


ModelSignal = Im1 + Im2 + Bkg;
