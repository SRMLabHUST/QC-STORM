
load('testdata.mat')
FitingData = double(YDat);

ParaNum=3;

% scaling factor, possibly can be automaticly determined
CoefA=10;
CoefS=1;

save ScalingCoef CoefA  CoefS

inInf = ParaPreEstimate(FitingData); % parameter preestimate

grad=NumericalGradient(inInf,FitingData); % get the gradient at current pameters
grad
D0=eye(ParaNum,ParaNum); % approximation of (Hessian)
I=eye(ParaNum,ParaNum);

circle=5; % iteration number

d0=-grad;

for i=1:circle 
    % normalize gradient
    td1=abs(d0);
    p1=mean(td1);
    
    d0=d0/p1*1.2;
        
    lpos=0.0001;
    rpos=1;
    
    ldat=LossFunction(inInf+d0*lpos, FitingData);
    rdat=LossFunction(inInf+d0*rpos, FitingData);
    
    
    for scnt=1:11 % find the best walk length by bisection method
        if(ldat<rdat)
            rpos=(lpos+rpos)/2;
            rdat=LossFunction(inInf+d0*rpos, FitingData);
        else
            lpos=(lpos+rpos)/2;
            ldat=LossFunction(inInf+d0*lpos, FitingData);
            
        end
    end
    pos=(lpos+rpos)/2;

    inInf=inInf+d0*pos;
    
    for ccnt=1:length(inInf)
        if(inInf(ccnt)<0)
            inInf(ccnt)=0; %parameter limitation
        end
    end        
%%
    sk=d0*pos;
    tgrad=grad;
    grad=NumericalGradient(inInf,FitingData);
    yk=grad-tgrad;
    % by BFGS method
    D0=(I-(sk*yk')/(yk'*sk))*D0*(I-(yk*sk')/(yk'*sk))+(sk*sk')/(yk'*sk);
    d0=-D0*grad;    %search direction

end



oInf=inInf
