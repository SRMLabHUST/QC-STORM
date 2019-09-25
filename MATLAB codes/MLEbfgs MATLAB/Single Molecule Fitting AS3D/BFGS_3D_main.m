% imput data, note target molecule should be (roughly) centered in the ROI
load('test_data_3D.mat')

InputROI = roi2;


% camera parameters
Offset = 398.6;
KAdc = 0.05;
QE = 0.9;


%%
[InInf] = BFGS_3D_f(InputROI, Offset, KAdc, QE);

%%
close all


InInf

figure
imshow(InputROI,[])
hold on
plot(InInf(2)+0.5, InInf(3)+0.5, 'bx','LineWidth',2, 'MarkerSize',8)

% figure
% imshow(RecoveredSignal, [])