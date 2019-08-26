function [fitresult, gof] = createFit_Sigma_ZDepth(SigmaDiff, ZDepth)
%CREATEFIT(SIGMADIFF1,ZDEPTH21)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : SigmaDiff1
%      Y Output: ZDepth21
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 03-Dec-2018 16:16:33 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( SigmaDiff, ZDepth );

% Set up fittype and options.
ft = fittype( 'poly4' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'ZDepth21 vs. SigmaDiff1', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel SigmaDiff1
% ylabel ZDepth21
% grid on
% 

