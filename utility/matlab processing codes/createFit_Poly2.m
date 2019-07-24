function [fitresult, gof] = createFit_Poly2(ZDepth2, sum2)
%CREATEFIT(ZDEPTH2,SUM2)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : ZDepth2
%      Y Output: sum2
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 26-Dec-2018 14:19:40 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( ZDepth2, sum2 );

% Set up fittype and options.
ft = fittype( 'poly2' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
% 
% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'sum2 vs. ZDepth2', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel ZDepth2
% ylabel sum2
% grid on


