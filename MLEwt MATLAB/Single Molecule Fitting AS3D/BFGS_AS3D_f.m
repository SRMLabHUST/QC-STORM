function [InInf] = BFGS_AS3D_f(InputROI, Offset, KAdc, QE, WLE_Enable)


%%
MoleculeSub = (double(InputROI) - Offset)*KAdc;


% moleculeSub = stdimg;

% scaling facor, large scaling factor for parameter with large danamic
% range, normolize the dynamic range to nearly [-1 1], large scaling factor enlarge the derivative
% if the fitting of some parameters converge slowly, probably you shold enlarge the scaling factor
ScalingCoef = [];
ScalingCoef.CoefA = 128;
ScalingCoef.CoefB = 32;
ScalingCoef.CoefS = 0.1;

ParaNum = 6; % fixed parameter number for 3D fitting

%%

[SigmaL, SigmaR, SigmaDiffH, SigmaU, SigmaD, SigmaDiffV] = GetWLEParaEstimation(MoleculeSub);

PSFSigmaX = min(SigmaL, SigmaR);
PSFSigmaY = min(SigmaU, SigmaD);

PSFWidth_Torlation	= 0.40;

% simple molecule clasification
if((SigmaDiffH > PSFWidth_Torlation)||(SigmaDiffV > PSFWidth_Torlation))
    PSFSigmaX = PSFSigmaX / 1.2;
    PSFSigmaY = PSFSigmaY / 1.2; 
else
    PSFSigmaX = 1000;
    PSFSigmaY = 1000;
end



WLEPara = [PSFSigmaX, PSFSigmaY];

%%
InInf = ParaPreEstimate_s3D(MoleculeSub, ScalingCoef); % parameter preestimate

grad = GetGradientNumerical_AS3D(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara); % get the gradient at current pameters

D0 = eye(ParaNum,ParaNum); % approximation of (Hessian)
I = eye(ParaNum,ParaNum);

circle = 11; % iteration number

d0 = -grad;

for i = 1:circle 


    td1 = abs(d0);
    p1 = mean(td1);
    
    d0 = d0/p1*1.2;
        
    lpos = 0.0001;
    rpos = 1;
    
    % find the best walk length by bisection method
    ldat = LossFunction_AS3D(InInf+d0*lpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    rdat = LossFunction_AS3D(InInf+d0*rpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    
    
    for scnt = 1:11
        if(ldat<rdat)
            rpos = (lpos+rpos)/2;
            rdat = LossFunction_AS3D(InInf+d0*rpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
        else
            lpos = (lpos+rpos)/2;
            ldat = LossFunction_AS3D(InInf+d0*lpos, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
            
        end
    end
    pos = (lpos+rpos)/2;

    InInf = InInf+d0*pos;
    
    % parameter limitation, only for parameters that far larger than 0
    % this is not very necessary since parameters are far from boundary
   for ccnt = 1:length(InInf)
        if(InInf(ccnt)<0)
            InInf(ccnt) = 0;
        end
    end        
    % update search direction by BFGS method
    sk = d0*pos;
    tgrad = grad;
    grad = GetGradientNumerical_AS3D(InInf, MoleculeSub, ScalingCoef, WLE_Enable, WLEPara);
    yk = grad-tgrad;
    D0 = (I-(sk*yk')/(yk'*sk))*D0*(I-(yk*sk')/(yk'*sk))+(sk*sk')/(yk'*sk);
    d0 = -D0*grad; % search direction
end




%% display recovered signal
ROISize = size(MoleculeSub,1);
RecoveredSignal = GetEstimatedSignal_AS3D(InInf, ScalingCoef, ROISize);
figure
imshow(RecoveredSignal, [])



%% remove scaling
InInf(1) = InInf(1)*ScalingCoef.CoefA/QE; % peak photon
InInf(6) = InInf(6)*ScalingCoef.CoefB/QE; % background intensity

InInf(4) = InInf(4)*ScalingCoef.CoefS; % Gaussian PSF sigma
InInf(5) = InInf(5)*ScalingCoef.CoefS; % Gaussian PSF sigma
InInf(4) = sqrt(0.5 / InInf(4));
InInf(5) = sqrt(0.5 / InInf(5));



