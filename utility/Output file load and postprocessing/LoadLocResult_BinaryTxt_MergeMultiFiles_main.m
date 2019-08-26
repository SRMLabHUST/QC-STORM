% to read localization point position list file
% how to use:
% copy the result file name to set following FileName, set current Matlab directory 
% to result folder directory, then run

clear,clc

FileName={
    'MMStack_Pos-1.ome_result2D7_M.txt',
    'MMStack_Pos-1_1.ome_result2D7_M.txt',
    'MMStack_Pos-1_2.ome_result2D7_M.txt',
    'MMStack_Pos-1_3.ome_result2D7_M.txt',
    'MMStack_Pos-1_4.ome_result2D7_M.txt',
    'MMStack_Pos-1_5.ome_result2D7_M.txt',
    'MMStack_Pos-1_6.ome_result2D7_M.txt',
    'MMStack_Pos-1_7.ome_result2D7_M.txt',
    'MMStack_Pos-1_8.ome_result2D7_M.txt',
    'MMStack_Pos-1_9.ome_result2D7_M.txt',
    'MMStack_Pos-1_10.ome_result2D7_M.txt',
    'MMStack_Pos-1_11.ome_result2D7_M.txt',
    'MMStack_Pos-1_12.ome_result2D7_M.txt',
    };

FileNum = length(FileName);

ParaNum = 12; % for 2d, 3d
% parameters are 
% peak intensity (photon), x (pixel), y (pixel), z (nm), PSFSigmaX (pixel), PSFSigmaY (pixel),Total intensity (photon), background (photon), SNR (peak to background e-), CRLBx (nm), CRLBy (nm), frame

LocArry = [];

FrameOffset = 0; % 0 if the first frame id is 1, or 1 if the first frame id is 0

for fcnt=1:FileNum
    curFileName=FileName{fcnt}
    

    fid=fopen(curFileName,'rb');
    if(fid==-1)
        break;
    end
    
    loc=fread(fid,inf,'float');
    fclose(fid);

    
    Len=floor(length(loc)/ParaNum);
    CurLocArry=zeros(ParaNum,Len);

    CurLocArry(:)=loc(1:ParaNum*Len);
    CurLocArry=CurLocArry';


    pos=CurLocArry(:,2)~=0;
    CurLocArry=CurLocArry(pos,:);
    CurLocArry=sortrows(CurLocArry,ParaNum);
    
    % frame index modification
    CurLocArry(:,end)=CurLocArry(:,end)+FrameOffset;
    FrameOffset=CurLocArry(end,end);
    
    LocArry=[LocArry; CurLocArry];
end


% plot(LocArry(:,2)+0.5,LocArry(:,3)+0.5,'x');

% length(CurLocArry(:,1))
TotalFrame = LocArry(end,end);
savename = sprintf("LocArry_whole_%df.mat", TotalFrame);

save(savename, 'LocArry' ,'-v7.3') 

FileName = sprintf("LocArry_whole_%df.txt", TotalFrame);

WriteLocResult_BinaryTxt(LocArry, FileName)
