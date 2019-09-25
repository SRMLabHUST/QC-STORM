function ModelSignal = EstimatedSignal_s2D(InInf, ScalingCoef, ROISize)

x=(1:ROISize)-0.5;
y=(1:ROISize)-0.5;

[X,Y] = meshgrid(x, y);


A1 = InInf(1)*ScalingCoef.CoefA;
x1 = InInf(2);
y1 = InInf(3);

Sigma1 = InInf(4)*ScalingCoef.CoefS;
Bkg = InInf(5)*ScalingCoef.CoefB;

Im1 = A1*exp(-((X-x1).*(X-x1) + (Y-y1).*(Y-y1))*Sigma1);


ModelSignal = Im1 + Bkg;

        