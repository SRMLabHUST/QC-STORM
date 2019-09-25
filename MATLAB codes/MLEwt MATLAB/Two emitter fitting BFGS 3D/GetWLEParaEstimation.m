function WLEPara = GetWLEParaEstimation(MoleculeSub)

ROISize = size(MoleculeSub,1);

if(ROISize<7)
    MarginLen = 0;
else
    MarginLen = 1;

end

rsel = MarginLen+1:ROISize-MarginLen;


MeanX = mean(MoleculeSub(rsel, rsel), 1);
MeanY = mean(MoleculeSub(rsel, rsel), 2)';

MeanX = MeanX - min(MeanX);
MeanY = MeanY - min(MeanY);

% plot(MeanX)
% hold on
% plot(MeanY)
[SigmaL, SigmaR, SigmaDiffH] = GetSigmaEstimation_Interp(MeanX);
[SigmaU, SigmaD, SigmaDiffV] = GetSigmaEstimation_Interp(MeanY);

PSFSigmaX = min(SigmaL, SigmaR);
PSFSigmaY = min(SigmaU, SigmaD);

PSFWidth_Torlation	= 0.20;

% simple molecule clasification

if (SigmaDiffH > PSFWidth_Torlation)
    PSFSigmaX = PSFSigmaX / 1.2;
else
    PSFSigmaX = PSFSigmaX * 1.2;
end

if (SigmaDiffV > PSFWidth_Torlation)
    PSFSigmaY = PSFSigmaY / 1.2; 
else
    PSFSigmaY = PSFSigmaY * 1.2;
end


WLEPara = [PSFSigmaX, PSFSigmaY];

