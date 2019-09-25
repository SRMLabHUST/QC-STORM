function [Deriv2_Pixels2X, Deriv2_Pixels2Y] = GaussianPSF_deriv2(x0, y0, Sigma, ROISize)

dh=0.01;

FluoD0 = Gaussian1(x0, y0, Sigma, ROISize);

FluoD0_dh_px = Gaussian1(x0+dh, y0, Sigma, ROISize);
FluoD0_dh_nx = Gaussian1(x0-dh, y0, Sigma, ROISize);

FluoD0_dh_py = Gaussian1(x0, y0+dh, Sigma, ROISize);
FluoD0_dh_ny = Gaussian1(x0, y0-dh, Sigma, ROISize);

% an array contain second derivation of each pixel to x and y
Deriv2_Pixels2X = (FluoD0_dh_px - 2*FluoD0 + FluoD0_dh_nx)/dh/dh;
Deriv2_Pixels2Y = (FluoD0_dh_py - 2*FluoD0 + FluoD0_dh_ny)/dh/dh;

