function [InInf] = BFGS_2D_f(InputROI, Offset, KAdc, QE)


%%
MoleculeSub = (double(InputROI) - Offset)*KAdc;


% the scaling coefficients to make the range of all parameters similar
ScalingCoef = [];

ScalingCoef.CoefA = 128;
ScalingCoef.CoefB = 32;
ScalingCoef.CoefS = 0.1;

ParaNum = 5; % fixed parameter number for 2D fitting


InInf = ParaPreEstimate_2D(MoleculeSub, ScalingCoef); % parameter preestimate

grad = NumericalGradient_2D(InInf, MoleculeSub, ScalingCoef); % get the gradient at current pameters

D0 = eye(ParaNum,ParaNum); % approximation of (Hessian)
I = eye(ParaNum,ParaNum);

circle = 8; % iteration number

d0 = -grad;

for i = 1:circle 


    td1 = abs(d0);
    p1 = mean(td1);
    
    d0 = d0/p1*1.2;
        
    lpos = 0.0001;
    rpos = 1;
    
    ldat = LossFunction_2D(InInf+d0*lpos, MoleculeSub, ScalingCoef);
    rdat = LossFunction_2D(InInf+d0*rpos, MoleculeSub, ScalingCoef);
    
    
    for scnt = 1:11 % find the best walk length by bisection method
        if(ldat<rdat)
            rpos = (lpos+rpos)/2;
            rdat = LossFunction_2D(InInf+d0*rpos, MoleculeSub, ScalingCoef);
        else
            lpos = (lpos+rpos)/2;
            ldat = LossFunction_2D(InInf+d0*lpos, MoleculeSub, ScalingCoef);
            
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
    grad = NumericalGradient_2D(InInf, MoleculeSub, ScalingCoef);
    yk = grad-tgrad;
    % by BFGS method
    D0 = (I-(sk*yk')/(yk'*sk))*D0*(I-(yk*sk')/(yk'*sk))+(sk*sk')/(yk'*sk);
    d0 = -D0*grad;    %search direction

end

%{
ROISize = size(MoleculeSub,1);
RecoveredSignal = EstimatedSignal_s2D(InInf, ScalingCoef, ROISize);
figure
imshow(RecoveredSignal, [])
%}

InInf(1) = InInf(1)*ScalingCoef.CoefA/QE; % peak photon
InInf(5) = InInf(5)*ScalingCoef.CoefB/QE; % background intensity

InInf(4) = InInf(4)*ScalingCoef.CoefS; % Gaussian PSF sigma
InInf(4) = sqrt(0.5 / InInf(4));
