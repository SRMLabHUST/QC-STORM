function [InInf, RecoveredSignal] = BFGS_3D_2E_f(InputROI, Offset, KAdc, QE, WLE_Enable, WLEPara)


%%
MoleculeSub = (double(InputROI) - Offset)*KAdc;


% the scaling coefficients to make the range of all parameters similar
ScalingCoef = [];

ScalingCoef.CoefA = 128;
ScalingCoef.CoefB = 32;
ScalingCoef.CoefS = 0.1;

ParaNum = 11;


ROISize = size(MoleculeSub,1);

InInf = ParaPreEstimate_3D_2E(MoleculeSub, ScalingCoef); % parameter preestimate


% WLEPara = GetWLEParaEstimation(MoleculeSub);

grad = possionfGradient_3D_2E(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara); % get the gradient at current pameters

D0 = eye(ParaNum,ParaNum); % approximation of (Hessian)
I = eye(ParaNum,ParaNum);

circle = 16; % iteration number

d0 = -grad;

for i = 1:circle 

    td1 = abs(d0);
    p1 = mean(td1);
    
    d0 = d0/p1*1.2;
        
    lpos = 0.0001;
    rpos = 1;
    
    ldat = possionf_3D_2E(InInf+d0*lpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    rdat = possionf_3D_2E(InInf+d0*rpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    
    
    for scnt = 1:11 % find the best walk length by bisection method
        if(ldat<rdat)
            rpos = (lpos+rpos)/2;
            rdat = possionf_3D_2E(InInf+d0*rpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
        else
            lpos = (lpos+rpos)/2;
            ldat = possionf_3D_2E(InInf+d0*lpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
            
        end
    end
    pos = (lpos+rpos)/2;

    InInf = InInf+d0*pos;
    
    for ccnt = 1:length(InInf)
        if(InInf(ccnt)<0)
            InInf(ccnt) = 0; %parameter limitation
        end
    end        
%%
    sk = d0*pos;
    tgrad = grad;
    grad = possionfGradient_3D_2E(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    yk = grad-tgrad;
    % by BFGS method
    D0 = (I-(sk*yk')/(yk'*sk))*D0*(I-(yk*sk')/(yk'*sk))+(sk*sk')/(yk'*sk);
    d0 = -D0*grad;    %search direction

end
InInf
RecoveredSignal = EstimatedSignal_3D_2E(InInf, ScalingCoef, ROISize);

figure
imshow(RecoveredSignal)

loss = sum(sum((MoleculeSub - RecoveredSignal).^2));

InInf([1 6]) = InInf([1 6])*ScalingCoef.CoefA/QE; % peak photon

InInf(11) = InInf(11)*ScalingCoef.CoefB/QE; % background intensity

InInf([4 5 9 10]) = InInf([4 5 9 10])*ScalingCoef.CoefS; % Gaussian PSF sigma
InInf([4 5 9 10]) = sqrt(0.5 ./ InInf([4 5 9 10]));
    


