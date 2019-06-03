function [SigmaL, SigmaR, SigmaDiffH, SigmaU, SigmaD, SigmaDiffV] = GetWLEParaEstimation(MoleculeSub)

ROISize = size(MoleculeSub,1);

MarginLen = 0;


rsel = MarginLen+1:ROISize-MarginLen;


MeanX = mean(MoleculeSub(rsel, rsel), 1);
MeanY = mean(MoleculeSub(rsel, rsel), 2)';

MeanX = MeanX - min(MeanX);
MeanY = MeanY - min(MeanY);

%{
plot(MeanX)
hold on
plot(MeanY)
%}
[SigmaL, SigmaR, SigmaDiffH] = GetSigmaEstimation_Interp(MeanX);
[SigmaU, SigmaD, SigmaDiffV] = GetSigmaEstimation_Interp(MeanY);


