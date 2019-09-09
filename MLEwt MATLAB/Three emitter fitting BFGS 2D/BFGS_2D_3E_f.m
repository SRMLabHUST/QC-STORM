function [InInf, RecoveredSignal] = BFGS_2D_3E_f(InputROI, Offset, KAdc, QE, WLE_Enable, WLEPara)


%%
MoleculeSub = (double(InputROI) - Offset)*KAdc;


% the scaling coefficients to make the range of all parameters similar
ScalingCoef = [];

ScalingCoef.CoefA = 128;
ScalingCoef.CoefB = 32;
ScalingCoef.CoefS = 0.1;

ParaNum = 11;


ROISize = size(MoleculeSub,1);

InInf = ParaPreEstimate_2D_3E(MoleculeSub, ScalingCoef); % parameter preestimate


% WLEPara = GetWLEParaEstimation(MoleculeSub);

grad = possionfGradient_2D_3E(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara); % get the gradient at current pameters

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
    
    ldat = possionf_2D_3E(InInf+d0*lpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    rdat = possionf_2D_3E(InInf+d0*rpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    
    
    for scnt = 1:11 % find the best walk length by bisection method
        if(ldat<rdat)
            rpos = (lpos+rpos)/2;
            rdat = possionf_2D_3E(InInf+d0*rpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
        else
            lpos = (lpos+rpos)/2;
            ldat = possionf_2D_3E(InInf+d0*lpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
            
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
    grad = possionfGradient_2D_3E(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    yk = grad-tgrad;
    % by BFGS method
    D0 = (I-(sk*yk')/(yk'*sk))*D0*(I-(yk*sk')/(yk'*sk))+(sk*sk')/(yk'*sk);
    d0 = -D0*grad;    %search direction

end

RecoveredSignal = EstimatedSignal_2D_3E(InInf, ScalingCoef, ROISize);

loss = sum(sum((MoleculeSub - RecoveredSignal).^2))

InInf([1 4 7]) = InInf([1 4 7])*ScalingCoef.CoefA/QE; % peak photon

InInf(11) = InInf(11)*ScalingCoef.CoefB/QE; % background intensity

InInf(10) = InInf(10)*ScalingCoef.CoefS; % Gaussian PSF sigma
InInf(10) = sqrt(0.5 / InInf(10));
    

