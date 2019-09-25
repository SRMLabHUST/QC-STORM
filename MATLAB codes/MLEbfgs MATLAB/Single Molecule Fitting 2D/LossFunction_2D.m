function VPoss = LossFunction_2D(InInf, MoleculeSub, ScalingCoef)


ROISize = size(MoleculeSub,1);

for i=1:length(InInf)
    if(InInf(i)<0)
        InInf(i)=0; % parameter restriction, only for parameter should larger than 0
    end
end

ModelSignal = EstimatedSignal_s2D(InInf, ScalingCoef, ROISize);


Loss = ModelSignal - MoleculeSub - MoleculeSub.*(log(ModelSignal) - log(max(MoleculeSub, 1)));


VPoss = sum(Loss(:));

