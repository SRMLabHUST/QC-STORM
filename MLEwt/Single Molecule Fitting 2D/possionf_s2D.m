function VPoss = possionf_s2D(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara)


ROISize = size(MoleculeSub,1);

for i=1:length(InInf)
    if(InInf(i)<0)
        InInf(i)=0;
    end
end

ModelSignal = EstimatedSignal_s2D(InInf, ScalingCoef, ROISize);


Loss = ModelSignal - MoleculeSub - MoleculeSub.*(log(ModelSignal) - log(max(MoleculeSub, 1)));

if(WLE_Enable)
    
    PSFSigmaX = WLEPara(1);
    PSFSigmaY = WLEPara(1);
    PSFSigmaX = min(PSFSigmaX, PSFSigmaY);
    PSFSigmaY = PSFSigmaX;
    
	WLE_Weitht = GetWeightMask(PSFSigmaX, PSFSigmaY, ROISize);

    Loss = Loss.*WLE_Weitht;
end

VPoss = sum(Loss(:));

