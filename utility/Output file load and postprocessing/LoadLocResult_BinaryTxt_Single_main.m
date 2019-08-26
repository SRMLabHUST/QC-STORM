% to read localization point position list file
% how to use:
% copy the result file name to set following FileName, set current Matlab directory 
% to result folder directory, then run

FileName = 'MMStack_Pos-1.ome-1_result3D13_S.txt';
ParaNum = 12; 

% parameters are 
% peak intensity (photon), x (pixel), y (pixel), z (nm), PSFSigmaX (pixel), PSFSigmaY (pixel), Total intensity (photon), background (photon), SNR (peak to background e-), CRLBx (nm), CRLBy (nm), frame

fid=fopen(FileName,'rb');
loc=fread(fid,inf,'float');
fclose(fid);

Len=floor(length(loc)/ParaNum);
LocArry=zeros(ParaNum,Len);

LocArry(:)=loc(1:ParaNum*Len);
LocArry=LocArry';

pos=LocArry(:,1)~=0;
LocArry=LocArry(pos,:);

LocArry=sortrows(LocArry,ParaNum);
% plot(LocArry(:,2)+0.5,LocArry(:,3)+0.5,'x');
% hold on

savename = sprintf("%s.mat", FileName(1:end-4));


save(savename, 'LocArry' ,'-v7.3') 
