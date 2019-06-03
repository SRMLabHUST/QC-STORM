function ModelSignal = EstimatedSignal_s3D(InInf, ScalingCoef, ROISize)

x=(1:ROISize)-0.5;
y=(1:ROISize)-0.5;

[X,Y] = meshgrid(x, y);


A1 = InInf(1)*ScalingCoef.CoefA;
x1 = InInf(2);
y1 = InInf(3);

SigmaX1 = InInf(4)*ScalingCoef.CoefS;
SigmaY1 = InInf(5)*ScalingCoef.CoefS;
Bkg = InInf(6)*ScalingCoef.CoefB;

Img1 = A1*exp(-((X-x1).*(X-x1)*SigmaX1 + (Y-y1).*(Y-y1)*SigmaY1));


ModelSignal = Img1 + Bkg;

        