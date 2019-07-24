close all


StartFrame = LocArry(1,end);
EndFrame = LocArry(end,end);
TotalFrame = EndFrame-StartFrame+1;

FluoPerFrame = zeros(1, TotalFrame);

for fcnt = 1:TotalFrame
    CurFrame = StartFrame + fcnt-1;
    
    pos1 = LocArry(:,end) == CurFrame;
    FluoPerFrame(fcnt) = sum(pos1);
end


plot(FluoPerFrame)

FrameSel1 = 13;
FrameSel2 = 30;

pos = (LocArry(:,end)>=FrameSel1)&(LocArry(:,end)<=FrameSel2);
LocArry = LocArry(pos,:);

figure
plot(LocArry(:,2), LocArry(:,3),'x');

save FilteredLocArray LocArry
