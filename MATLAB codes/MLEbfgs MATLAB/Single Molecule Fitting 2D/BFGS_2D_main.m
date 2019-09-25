
% imput data, note target molecule should be (roughly) centered in the ROI
load('testdata.mat')

InputROI = roi1;


% camera parameters
Offset = 100;
KAdc = 0.45;
QE = 0.72;


%%


[InInf] = BFGS_2D_f(InputROI, Offset, KAdc, QE);



%%
close all

InInf

figure
imshow(InputROI,[])
hold on
plot(InInf(2)+0.5, InInf(3)+0.5, 'bx','LineWidth',2, 'MarkerSize',8)


