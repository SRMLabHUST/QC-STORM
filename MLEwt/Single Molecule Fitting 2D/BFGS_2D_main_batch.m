
% imput data, note target molecule should be (roughly) centered in the ROI


load('MoleculeDistance_Simu.mat')

FrameNum = 8000;

% find the center position
% image1 = mean(Image_All,3);
% imwrite(uint16(image1), 'image1.tif')

ROIPosX = 16 +1;
ROIPosY = 24 +1;

ROISize = 11;
ROISize_Half = floor(ROISize/2);

% camera parameters
Offset = 100;
KAdc = 0.45;
QE = 0.72;

% choose between MLE and WLE
WLE_Enable = 0;

CenterPosX = ROISize/2;
CenterPosY = ROISize/2;

DatLen = length(MoleculeDistance_nm);

RMSE_Vary = zeros(1,DatLen);

for i = 1:DatLen
    disp(i)
    
    CurDistance = MoleculeDistance_nm(i);
    
    InInf_all = zeros(FrameNum,5);
    
   	savename = sprintf('Image_Distance%d.mat', CurDistance);
    load(savename)
    
    for fcnt = 1:FrameNum
%         ImageName = sprintf('ImageSimu_MoleculePair_Distance%d/img%d.tif', CurDistance, fcnt);
%         
%         curImage = imread(ImageName);
        
        curImage = Image_All(:, :, fcnt);

        rsel = ROIPosY-ROISize_Half : ROIPosY+ROISize_Half;
        csel = ROIPosX-ROISize_Half : ROIPosX+ROISize_Half;
        
        InputROI = curImage(rsel, csel);
        
%         imshow(InputROI,[]);
        
        InInf_t = BFGS_2D_f(InputROI, Offset, KAdc, QE, WLE_Enable);
    
        InInf_all(fcnt,:) = InInf_t';
    end
    
    
    diff = (InInf_all(:,2)-CenterPosX).^2 + (InInf_all(:,3)-CenterPosY).^2;
    RMSE = sqrt(mean(diff));
    
%     savename = sprintf('Ininf_Distance%d_Single_WLE%d.mat', CurDistance, WLE_Enable);
%     save(savename,'InInf_all', 'RMSE', 'ROIPosX', 'ROIPosY', 'CenterPosX', 'CenterPosY')
% 
    
    RMSE_Vary(i) = RMSE;
end


savename = sprintf('RMSE_Vary_Single_WLE%d_ROISize%d.mat', WLE_Enable, ROISize);

save(savename, 'RMSE_Vary','MoleculeDistance_nm')  

plot(MoleculeDistance_nm, RMSE_Vary, '-x', 'LineWidth', 2, 'MarkerSize', 8)
hold on

hold on

xlabel('Distance (nm)','FontSize', 16);
ylabel('RMSE (nm)', 'FontSize', 16); % 

ha = gca; % current axis
set(ha,'FontSize', 16, 'LineWidth', 2);

legend('GUEST off 7x7', 'GUEST off 11x11', 'GUEST on 7x7', 'GUEST on 11x11')

