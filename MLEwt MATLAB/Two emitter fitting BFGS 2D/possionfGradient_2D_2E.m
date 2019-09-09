function grad = possionfGradient_2D_2E(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara)

for i  =  1:length(InInf)
    if(InInf(i)<0)
        InInf(i) = 0; %
    end
end

nx = length(InInf);
grad = zeros(nx,1);

for i = 1:nx
	b1 = InInf;
	b2 = InInf;
    
	b1(i) = InInf(i)-0.001;
    L1 = possionf_2D_2E(b1, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    
	b2(i) = InInf(i)+0.001;
    L2 = possionf_2D_2E(b2, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
	
    if(InInf(i)<0)
       grad(i) = 0; 
    else
        grad(i) = (L2-L1)/0.002;
    end
end
