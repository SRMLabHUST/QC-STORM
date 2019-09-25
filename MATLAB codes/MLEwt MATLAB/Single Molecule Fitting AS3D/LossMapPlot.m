load('test_data_3D.mat')

InputROI = roi0;

WLEEn = 1;

% camera parameters
Offset = 398.6;
KAdc = 0.05;
QE = 0.9;

MoleculeSub = (double(InputROI) - Offset)*KAdc;

ROISize = size(InputROI,1);

InInf = [
  170.1156
    5.8962
    5.3566
    1.2551
    1.7071
   39.1890
   ];

ScalingCoef.CoefA = 1;
ScalingCoef.CoefB = 1;
ScalingCoef.CoefS = 1;

x=0:0.1:ROISize;
y=0:0.1:ROISize;
datlen = length(x);

LossData=zeros(datlen,datlen);

for ycnt=1:datlen
    for xcnt=1:datlen
        
        InInf(2) = x(xcnt);
        InInf(3) = y(ycnt);
        
        ModelSignal = EstimatedSignal_s3D(InInf, ScalingCoef, ROISize);


        Loss = ModelSignal - MoleculeSub - MoleculeSub.*(log(ModelSignal) - log(max(MoleculeSub, 1)));
        
        if(WLEEn)
            [SigmaL, SigmaR, SigmaDiffH, SigmaU, SigmaD, SigmaDiffV] = GetWLEParaEstimation(MoleculeSub);

            PSFSigmaX = min(SigmaL, SigmaR);
            PSFSigmaY = min(SigmaU, SigmaD);

            PSFWidth_Torlation	= 0.40;

            % simple molecule clasification
            if((SigmaDiffH > PSFWidth_Torlation)||(SigmaDiffV > PSFWidth_Torlation))
                PSFSigmaX = PSFSigmaX / 1.2;
                PSFSigmaY = PSFSigmaY / 1.2; 
            else
                PSFSigmaX = 1000;
                PSFSigmaY = 1000;
            end

            WLE_Weitht = GetWeightMask(PSFSigmaX, PSFSigmaY, ROISize);

            Loss = Loss.*WLE_Weitht;
        end
        VPoss = sum(Loss(:));      
        
        LossData(ycnt,xcnt)=VPoss;
    end
end
figure
imshow(LossData,[])
% mesh(x, y, LossData)

 xlabel('X','FontSize', 16);
 ylabel('Y', 'FontSize', 16);

ha = gca; % current axis
set(ha, 'FontSize', 16, 'LineWidth', 2);
