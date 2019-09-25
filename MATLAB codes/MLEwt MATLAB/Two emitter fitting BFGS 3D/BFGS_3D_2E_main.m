
% imput data, note target molecule should be (roughly) centered in the ROI
% load('test_data_2D.mat')

InputROI = roi2;


% InputROI = rot90(InputROI);

ROISize = size(InputROI,1);

% camera parameters
Offset = 100;
KAdc = 0.45;
QE = 0.9;

% choose between MLE and WLE
WLE_Enable = 0;

WLEPara = [ROISize/1.9/2.35, ROISize/1.9/2.35];

%%
InInf = BFGS_3D_2E_f(InputROI, Offset, KAdc, QE, WLE_Enable, WLEPara);

%%

InInf

figure
imshow(InputROI,[])
hold on
plot(InInf(2), InInf(3), 'bx','LineWidth',2, 'MarkerSize',8)
plot(InInf(7), InInf(8), 'bx','LineWidth',2, 'MarkerSize',8)

