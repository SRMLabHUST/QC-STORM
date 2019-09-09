function InInf = BFGS_2D_2E_f(InputROI, Offset, KAdc, QE, WLE_Enable, WLEPara)


%%
MoleculeSub = (double(InputROI) - Offset)*KAdc;


% the scaling coefficients to make the range of all parameters similar
ScalingCoef = [];

ScalingCoef.CoefA = 128;
ScalingCoef.CoefB = 64;
ScalingCoef.CoefS = 0.1;

ParaNum = 8;


ROISize = size(MoleculeSub,1);

InInf = ParaPreEstimate_2D_2E(MoleculeSub, ScalingCoef); % parameter preestimate

% InInf(2)= 4;
% InInf(3)= 4;
% 
% InInf(5)= 6;
% InInf(6)= 4;


% WLEPara = GetWLEParaEstimation(MoleculeSub);

grad = possionfGradient_2D_2E(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara); % get the gradient at current pameters

D0 = eye(ParaNum,ParaNum); % approximation of (Hessian)
I = eye(ParaNum,ParaNum);

circle = 14; % iteration number

d0 = -grad;

for i = 1:circle 


    td1 = abs(d0);
    p1 = mean(td1);
    
    d0 = d0/p1*1.2;
        
    lpos = 0.0001;
    rpos = 1;
    
    ldat = possionf_2D_2E(InInf+d0*lpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    rdat = possionf_2D_2E(InInf+d0*rpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    
    
    for scnt = 1:11 % find the best walk length by bisection method
        if(ldat<rdat)
            rpos = (lpos+rpos)/2;
            rdat = possionf_2D_2E(InInf+d0*rpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
        else
            lpos = (lpos+rpos)/2;
            ldat = possionf_2D_2E(InInf+d0*lpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
            
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
    grad = possionfGradient_2D_2E(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    yk = grad-tgrad;
    % by BFGS method
    D0 = (I-(sk*yk')/(yk'*sk))*D0*(I-(yk*sk')/(yk'*sk))+(sk*sk')/(yk'*sk);
    d0 = -D0*grad;    %search direction

end

% RecoveredSignal = EstimatedSignal_2D_2E(InInf, ScalingCoef, ROISize);

InInf(1) = InInf(1)*ScalingCoef.CoefA/QE; % peak photon
InInf(4) = InInf(4)*ScalingCoef.CoefA/QE; % peak photon

InInf(8) = InInf(8)*ScalingCoef.CoefB/QE; % background intensity

InInf(7) = InInf(7)*ScalingCoef.CoefS; % Gaussian PSF sigma
InInf(7) = sqrt(0.5 / InInf(7));



