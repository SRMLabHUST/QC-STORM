function VPoss = LossFunction_AS3D(InInf, MoleculeSub, ScalingCoef)

for ccnt = 1:length(InInf)
    if(InInf(ccnt)<0)
        InInf(ccnt) = 0; % parameter restriction, only for parameter should larger than 0
    end
end        

ROISize = size(MoleculeSub,1);

ModelSignal = EstimatedSignal_s3D(InInf, ScalingCoef, ROISize);


Loss = ModelSignal - MoleculeSub - MoleculeSub.*(log(ModelSignal) - log(max(MoleculeSub, 1)));


VPoss = sum(Loss(:));



