function grad=NumericalGradient(inInf, FitingData)


for i=1:length(inInf)
    if(inInf(i)<0)
        inInf(i)=0; %范围限制
    end
end

nx=length(inInf);
grad=zeros(nx,1);

for i=1:nx
    %导数的近似方法求各个分量的偏导
	b1=inInf;
	b2=inInf;
    
	b1(i)=inInf(i)-0.001;
    L1=LossFunction(b1, FitingData);
    
	b2(i)=inInf(i)+0.001;
    L2=LossFunction(b2, FitingData);
    
    %不能用单边微分，否则计算不够准确，已经实际验证
	
    if(inInf(i)<0)
       grad(i)=0; 
    else
        grad(i)=(L2-L1)/0.002;
    end
end
