
% imput data, note target molecule should be (roughly) centered in the ROI
load('test_data_2D.mat')

InputROI = roi3;


% camera parameters
Offset = 100;
KAdc = 0.45;
QE = 0.72;

% choose between MLE and WLE
WLE_Enable = 1;


%%
[InInf] = BFGS_2D_f(InputROI, Offset, KAdc, QE, WLE_Enable);

%%

InInf

figure
imshow(InputROI,[])
hold on
plot(InInf(2)+0.5, InInf(3)+0.5, 'bx','LineWidth',2, 'MarkerSize',8)

