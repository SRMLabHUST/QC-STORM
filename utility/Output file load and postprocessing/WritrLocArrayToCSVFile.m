
PixelSizeX = 117; % nm
PixelSizeY = 127; % nm

FileName = 'Tubulin-A647-3D beads.csv';

LocArry(:,2) = LocArry(:,2)*PixelSizeX;
LocArry(:,3) = LocArry(:,3)*PixelSizeY;

% peak intensity (photon), x (pixel), y (pixel), z (nm), PSFSigmaX (pixel), PSFSigmaY (pixel),Total intensity (photon), background (photon), SNR (peak to background e-), CRLBx (nm), CRLBy (nm), frame

ColName = {'PeakPhoton', 'Xnano', 'Ynano', 'Znano', 'PSFSigmaX', 'PSFSigmaY', 'TotalPhoton','BackgroundPhoton', 'pSNR', 'CRLBX', 'CRLBY','Frame'};

PeakPhoton = LocArry(:,1);
XNano = LocArry(:,2);
YNano = LocArry(:,3);
ZNano = LocArry(:,4);
PSFSigmaX = LocArry(:,5);
PSFSigmaY = LocArry(:,6);
TotalPhoton = LocArry(:,7);
Background = LocArry(:,8);
pSNR = LocArry(:,9);
CRLBX = LocArry(:,10);
CRLBY = LocArry(:,11);
Frame = LocArry(:,12);



data = table(PeakPhoton, XNano, YNano, ZNano, PSFSigmaX, PSFSigmaY, TotalPhoton, Background, pSNR, CRLBX, CRLBY, Frame, 'VariableNames', ColName); 

writetable(data, FileName)
