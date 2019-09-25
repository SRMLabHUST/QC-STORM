function oinf= ParaPreEstimate(FitingData)

% CoefA=256;
% CoefB=128;
% CoefS=0.1;
load ScalingCoef

Amp=98;
X0=38;
Sigma=6;

oinf = [Amp/CoefA X0 Sigma/CoefS]';


% RegionLength=size(moleculeSub,1);
% 
% 
% imgdat=double(moleculeSub);
% backg=mean([imgdat(1,1:end) imgdat(end,1:end) imgdat(2:end-1,1)' imgdat(2:end-1,end)']);
% 
% Amplitude=max(imgdat(:))-backg;
% 
% X0=RegionLength/2;
% Y0=RegionLength/2;
% SigmaX = 1.4;
% 
% sigmax=1/(2*SigmaX*SigmaX);
% 
% Amplitude=Amplitude/CoefA;
% backg=backg/CoefB;
% sigmax=sigmax/CoefS;
% 
% oinf=[Amplitude; X0; Y0; sigmax; backg];
% 

