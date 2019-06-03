function WLE_Weitht = GetWeightMask(PSFSigmaX, PSFSigmaY, ROISize)

    

CenterPos = ROISize/2;

x = (1:ROISize)-0.5;
y = (1:ROISize)-0.5;

[X,Y] = meshgrid(x, y);
WLE_Weitht = exp(-((X-CenterPos).*(X-CenterPos)/(2*PSFSigmaX*PSFSigmaX) + (Y-CenterPos).*(Y-CenterPos)/(2*PSFSigmaY*PSFSigmaY)));

% mesh(WLE_Weitht)
