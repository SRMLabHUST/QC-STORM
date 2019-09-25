function VPoss = LossFunction(inInf, FitingData)

% keep the same with outside
% CoefA=256;
% CoefB=128;
% CoefS=0.1;
load ScalingCoef

for i=1:length(inInf)
    if(inInf(i)<0)
        inInf(i)=0; %·¶Î§ÏÞÖÆ
    end
end

Amp=inInf(1);
X0=inInf(2);
Sigma=inInf(3);


DatLen=length(FitingData);

Sigma1=1/(2*Sigma*CoefS*Sigma*CoefS);

XPos = (1:1:DatLen) - 0.5;
YDat = Amp*CoefA * exp(-(XPos - X0).^2 * Sigma1);

VPoss=sum((YDat-FitingData).^2);




