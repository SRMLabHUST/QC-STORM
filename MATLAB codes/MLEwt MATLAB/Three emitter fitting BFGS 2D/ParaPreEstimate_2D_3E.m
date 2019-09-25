function oInf = ParaPreEstimate_2D_3E(InputROI, ScalingCoef)


ROISize = size(InputROI,1);
ROISize_H = round(ROISize/2);


imgdat = double(InputROI);
Bkg = mean([imgdat(1,1:end) imgdat(end,1:end) imgdat(2:end-1,1)' imgdat(2:end-1,end)']);

Amp = max(imgdat(:)) - Bkg;

SigmaX = ROISize/2/2.35;
SigmaX1 = 1/(2*SigmaX*SigmaX);

%% scaling
Amp = Amp/ScalingCoef.CoefA;
Bkg = Bkg/ScalingCoef.CoefB;

SigmaX1 = SigmaX1/ScalingCoef.CoefS;

%%
xpos_ini = zeros(1,3);
ypos_ini = zeros(1,3);

SumData = zeros(1,4);

SumData(1) = sum(sum(InputROI(1:ROISize_H, 1:ROISize_H)));
SumData(2) = sum(sum(InputROI(1:ROISize_H, ROISize_H:end)));
SumData(3) = sum(sum(InputROI(ROISize_H:end, 1:ROISize_H)));
SumData(4) = sum(sum(InputROI(ROISize_H:end, ROISize_H:end)));

mindat = min(SumData);

gap = ROISize/4.6;

ValidPos = 1;
for i=1:4
    if(SumData(i)~=mindat)
        xpos_ini(ValidPos) = ROISize/2 + (-1)^(mod(i,2))*gap;
        ypos_ini(ValidPos) = ROISize/2 + (-1)^(floor((i+1)/2))*gap;
        
        ValidPos = ValidPos+1;
    end
end

% figure
% imshow(InputROI,[])
% hold on
% plot(xpos_ini,ypos_ini,'yx')

%%
A1 = 0.8*Amp;
A2 = 0.8*Amp;
A3 = 0.8*Amp;

x1 = xpos_ini(1);
y1 = ypos_ini(1);

x2 = xpos_ini(2);
y2 = ypos_ini(2);

x3 = xpos_ini(3);
y3 = ypos_ini(3);


oInf = [A1; x1; y1; A2; x2; y2; A3; x3; y3; SigmaX1; Bkg];


