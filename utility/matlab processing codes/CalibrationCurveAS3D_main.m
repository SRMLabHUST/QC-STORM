% load('FilteredLocArray.mat')

%% generate Calibration Data from localizations of beads stack

close all

FrameDepthGap = 10; % 2*15.625; % nm

FocalPlaneFrameChose_Manual = 0;
FocalPlaneFrame_Manual = 21;

FrameOrder = 0; % 0: negative depth to positive depth, 1: positive depth to negative depth



%%
FrameSel = LocArry(1,end):1:LocArry(end,end);
TotalFrame = length(FrameSel);


Distance_Th = 0.5;

SigmaXVary = zeros(1, TotalFrame);
SigmaYVary = zeros(1, TotalFrame);

for fcnt = 1:TotalFrame
    CurFrame = FrameSel(fcnt);
    
    pos = LocArry(:,end) == CurFrame;
    LocArry1 = LocArry(pos,:);
    
    
    SigmaXY_raw = LocArry1(:,[5 6]);
    DatLen = size(SigmaXY_raw,1);

    CenterPos = mean(SigmaXY_raw,1);
    for ItNum = 1:6
        diff = SigmaXY_raw - repmat(CenterPos, DatLen,1);
        distance = sqrt(sum(diff.^2,2));
        pos = distance<Distance_Th;
        SigmaXY=SigmaXY_raw(pos,:);
       	CenterPos = mean(SigmaXY,1);

    end
%     plot(SigmaXY_raw(:,1),SigmaXY_raw(:,2),'x')
%     hold on
%     plot(SigmaXY(:,1),SigmaXY(:,2),'x')
%     plot(CenterPos(:,1),CenterPos(:,2),'o')

    SigmaXVary(fcnt) = CenterPos(:,1);
    SigmaYVary(fcnt) = CenterPos(:,2);

end

diff = abs(SigmaXVary - SigmaYVary);
[~,MinPos] = min(diff);

FocalPlaneFrame = FrameSel(MinPos);

if(FocalPlaneFrameChose_Manual)
    FocalPlaneFrame = FocalPlaneFrame_Manual;
end



figure
plot(FrameSel, SigmaXVary, 'x', 'LineWidth',2, 'MarkerSize',8)
hold on
plot(FrameSel, SigmaYVary, 'x', 'LineWidth',2, 'MarkerSize',8)
legend('x','y')



ZDepth2 = (FrameSel - FocalPlaneFrame) * FrameDepthGap;

if(FrameOrder > 0)
    ZDepth2 = (FocalPlaneFrame - FrameSel) * FrameDepthGap;
end


SigmaDiff = SigmaXVary.^2 - SigmaYVary.^2;

figure
plot(ZDepth2, SigmaXVary, 'x', 'LineWidth',2, 'MarkerSize',8)
hold on
plot(ZDepth2, SigmaYVary, 'x', 'LineWidth',2, 'MarkerSize',8)
legend('x','y')

posL = SigmaXVary >= SigmaYVary - 0.2;
posR = SigmaXVary <= SigmaYVary + 0.2;


ZDepth2_L = ZDepth2(posL);
SigmaDiff_L = SigmaDiff(posL);

ZDepth2_R = ZDepth2(posR);
SigmaDiff_R = SigmaDiff(posR);



[fitresult_L, ~] = createFit_Sigma_ZDepth(SigmaDiff_L, ZDepth2_L);
[fitresult_R, ~] = createFit_Sigma_ZDepth(SigmaDiff_R, ZDepth2_R);


Calib_SigmaXGSigmaY.p4 = fitresult_L.p1;
Calib_SigmaXGSigmaY.p3 = fitresult_L.p2;
Calib_SigmaXGSigmaY.p2 = fitresult_L.p3;
Calib_SigmaXGSigmaY.p1 = fitresult_L.p4;
Calib_SigmaXGSigmaY.p0 = fitresult_L.p5;

Calib_SigmaXLSigmaY.p4 = fitresult_R.p1;
Calib_SigmaXLSigmaY.p3 = fitresult_R.p2;
Calib_SigmaXLSigmaY.p2 = fitresult_R.p3;
Calib_SigmaXLSigmaY.p1 = fitresult_R.p4;
Calib_SigmaXLSigmaY.p0 = fitresult_R.p5;

Calib_SigmaXGSigmaY
Calib_SigmaXLSigmaY

figure

plot(SigmaDiff, ZDepth2, 'x', 'LineWidth',2, 'MarkerSize',8);
hold on
ZDepth2_f_L = Calib_SigmaXGSigmaY.p4*SigmaDiff_L.^4 + Calib_SigmaXGSigmaY.p3*SigmaDiff_L.^3 + Calib_SigmaXGSigmaY.p2*SigmaDiff_L.^2 + Calib_SigmaXGSigmaY.p1*SigmaDiff_L + Calib_SigmaXGSigmaY.p0;
plot(SigmaDiff_L, ZDepth2_f_L, '-o', 'LineWidth',2, 'MarkerSize',8);

ZDepth2_f_R = Calib_SigmaXLSigmaY.p4*SigmaDiff_R.^4 + Calib_SigmaXLSigmaY.p3*SigmaDiff_R.^3 + Calib_SigmaXLSigmaY.p2*SigmaDiff_R.^2 + Calib_SigmaXLSigmaY.p1*SigmaDiff_R + Calib_SigmaXGSigmaY.p0;
plot(SigmaDiff_R, ZDepth2_f_R, '-o', 'LineWidth',2, 'MarkerSize',8);


xlabel('SigmaX^2 - SigmaY^2', 'FontSize', 16);
ylabel('Z depth (nm)','FontSize', 16);
ha = gca; % current axis
set(ha, 'FontSize', 16, 'LineWidth', 2);


save CalibrationData3D fitresult_L Calib_SigmaXGSigmaY Calib_SigmaXLSigmaY fitresult_R SigmaDiff ZDepth2 FrameDepthGap FocalPlaneFrame FrameSel SigmaXVary SigmaYVary FrameOrder

