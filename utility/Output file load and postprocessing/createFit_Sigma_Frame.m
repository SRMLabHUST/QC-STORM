function [fitresult, gof] = createFit_Sigma_Frame(FrameId, SigmaYVary)
%CREATEFIT(FRAMEID,SIGMAYVARY)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : FrameId
%      Y Output: SigmaYVary
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 03-Dec-2018 16:00:49 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( FrameId, SigmaYVary );

% Set up fittype and options.
ft = fittype( 'poly4' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'SigmaYVary vs. FrameId', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel FrameId
% ylabel SigmaYVary
% grid on


