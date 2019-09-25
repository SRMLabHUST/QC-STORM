function [SigmaL, SigmaR, SigmaDiff] = GetSigmaEstimation_Interp(MeanData)


ROISize = length(MeanData);

MeanData = MeanData-min(MeanData);

x=1:ROISize;
% plot(x, MeanData)
% hold on

IntGap = 0.1;
x1=1:IntGap:ROISize;

MeanData1 = interp1(x, MeanData, x1,'spline');

% MeanData1 = SplineInterpolation(x, MeanData, x1);

%{
plot(x1, MeanData1,'-x')
hold on
plot(x, MeanData,'-o')
%}

Center_Interplot = round(length(MeanData1)/2);
curx=0;
for i = Center_Interplot-10:Center_Interplot+10
    if(MeanData1(i)>curx)
        curx = MeanData1(i);
        CenterPos_Int = i;
    end
end

Peak_Rough = MeanData1(CenterPos_Int);

LLen = CenterPos_Int-1;
RLen = length(MeanData1)-CenterPos_Int;
CmpLen = min(LLen,RLen);


%% left side

SDat_L = 0;
for i=(CenterPos_Int-CmpLen):(CenterPos_Int) %-1
    SDat_L = SDat_L + MeanData1(i);
end


%% right side
SDat_R = 0;
for i=(CenterPos_Int):(CenterPos_Int+CmpLen) %+1
    SDat_R = SDat_R + MeanData1(i);
end

%%


SigmaL = SDat_L*IntGap/Peak_Rough/sqrt(pi/2);
SigmaR = SDat_R*IntGap/Peak_Rough/sqrt(pi/2);

SigmaDiff = abs(SigmaL - SigmaR)/min(SigmaL, SigmaR);


