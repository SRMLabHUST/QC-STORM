function grad = NumericalGradient_2D(InInf, MoleculeSub, ScalingCoef)



nx = length(InInf);
grad = zeros(nx,1);

for i = 1:nx
	b1 = InInf;
	b2 = InInf;
    
	b1(i) = InInf(i)-0.001;
    L1 = LossFunction_2D(b1, MoleculeSub, ScalingCoef);
	b2(i) = InInf(i)+0.001;
    L2 = LossFunction_2D(b2, MoleculeSub, ScalingCoef);
	
    grad(i) = (L2-L1)/0.002;
end
