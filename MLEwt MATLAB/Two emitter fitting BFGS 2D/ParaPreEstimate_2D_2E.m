function oInf = ParaPreEstimate_2D_2E(InputROI, ScalingCoef)


ROISize = size(InputROI,1);
ROISize_H = round(ROISize/2);


imgdat=double(InputROI);
Bkg = mean([imgdat(1,1:end) imgdat(end,1:end) imgdat(2:end-1,1)' imgdat(2:end-1,end)']);

Amp = max(imgdat(:)) - Bkg;

SigmaX = ROISize/2/2.35;
SigmaX1 = 1/(2*SigmaX*SigmaX);

%% scaling
Amp = Amp/ScalingCoef.CoefA;
Bkg = Bkg/ScalingCoef.CoefB;

SigmaX1 = SigmaX1/ScalingCoef.CoefS;

%%

SumData = zeros(1,4);

SumData(1) = sum(sum(InputROI(1:ROISize_H, 1:ROISize_H)));
SumData(2) = sum(sum(InputROI(1:ROISize_H, ROISize_H:end)));
SumData(3) = sum(sum(InputROI(ROISize_H:end, 1:ROISize_H)));
SumData(4) = sum(sum(InputROI(ROISize_H:end, ROISize_H:end)));

[~,MaxPos] = max(SumData);

gap = ROISize/4.6;




%%
A1 = 0.8*Amp;
A2 = 0.8*Amp;

x1 = ROISize/2 ;
y1 = ROISize/2 ;
x2 = ROISize/2 + (-1)^(mod(MaxPos,2))*gap;
y2 = ROISize/2 + (-1)^(floor((MaxPos+1)/2))*gap;


oInf = [A1; x1; y1; A2; x2; y2; SigmaX1; Bkg];


