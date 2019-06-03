function VPoss = possionf_s3D(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara)


ROISize = size(MoleculeSub,1);

for i=1:length(InInf)
    if(InInf(i)<0)
        InInf(i)=0;
    end
end

ModelSignal = EstimatedSignal_s3D(InInf, ScalingCoef, ROISize);


Loss = ModelSignal - MoleculeSub - MoleculeSub.*(log(ModelSignal) - log(max(MoleculeSub, 1)));

if(WLE_Enable)

    PSFSigmaX = WLEPara(1);
    PSFSigmaY = WLEPara(2);
    
	WLE_Weitht = GetWeightMask(PSFSigmaX, PSFSigmaY, ROISize);
        
    Loss = Loss.*WLE_Weitht;
end

VPoss = sum(Loss(:));



