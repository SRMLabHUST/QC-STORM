function [InInf] = BFGS_2D_f(InputROI, Offset, KAdc, QE, WLE_Enable)


%%
MoleculeSub = (double(InputROI) - Offset)*KAdc;


% the scaling coefficients to make the range of all parameters similar
ScalingCoef = [];

ScalingCoef.CoefA = 128;
ScalingCoef.CoefB = 32;
ScalingCoef.CoefS = 0.1;

ParaNum = 5; % fixed parameter number for 2D fitting



%% WLE parameter estimation

[SigmaL, SigmaR, SigmaDiffH, SigmaU, SigmaD, SigmaDiffV] = GetWLEParaEstimation(MoleculeSub);


PSFSigmaX = min(SigmaL, SigmaR);
PSFSigmaY = min(SigmaU, SigmaD);

PSFWidth_Torlation	= 0.40;

% simple molecule clasification


if((SigmaDiffH > PSFWidth_Torlation)||(SigmaDiffV > PSFWidth_Torlation))
    % is contaminated
    PSFSigmaX = PSFSigmaX / 1.2;
    PSFSigmaY = PSFSigmaY / 1.2; 
else
    PSFSigmaX = 1000;
    PSFSigmaY = 1000;
end



WLEPara = [PSFSigmaX, PSFSigmaY];

%%

InInf = ParaPreEstimate_s2D(MoleculeSub, ScalingCoef); % parameter preestimate

grad = possionfGradient_s2D(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara); % get the gradient at current pameters

D0 = eye(ParaNum,ParaNum); % approximation of (Hessian)
I = eye(ParaNum,ParaNum);

circle = 9; % iteration number

d0 = -grad;

for i = 1:circle 


    td1 = abs(d0);
    p1 = mean(td1);
    
    d0 = d0/p1*1.2;
        
    lpos = 0.0001;
    rpos = 1;
    
    ldat = possionf_s2D(InInf+d0*lpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    rdat = possionf_s2D(InInf+d0*rpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    
    
    for scnt = 1:11 % find the best walk length by bisection method
        if(ldat<rdat)
            rpos = (lpos+rpos)/2;
            rdat = possionf_s2D(InInf+d0*rpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
        else
            lpos = (lpos+rpos)/2;
            ldat = possionf_s2D(InInf+d0*lpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
            
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
    grad = possionfGradient_s2D(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    yk = grad-tgrad;
    % by BFGS method
    D0 = (I-(sk*yk')/(yk'*sk))*D0*(I-(yk*sk')/(yk'*sk))+(sk*sk')/(yk'*sk);
    d0 = -D0*grad;    %search direction

end


ROISize = size(MoleculeSub,1);
RecoveredSignal = EstimatedSignal_s2D(InInf, ScalingCoef, ROISize);
figure
imshow(RecoveredSignal, [])

%% remove scaling

InInf(1) = InInf(1)*ScalingCoef.CoefA/QE; % peak photon
InInf(5) = InInf(5)*ScalingCoef.CoefB/QE; % background intensity

InInf(4) = InInf(4)*ScalingCoef.CoefS; % Gaussian PSF sigma
InInf(4) = sqrt(0.5 / InInf(4));
