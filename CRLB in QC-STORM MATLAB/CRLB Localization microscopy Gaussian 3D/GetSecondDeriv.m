
function [deri2, deri2_ln] = GetSecondDeriv(InputPara, ROISize)

% InputPara=[Amp,x0,y0,Sigma,Bkg]

% iA=100;
% iX0=0;
% iSigma=1;

dh = 0.02;

ParaNum = length(InputPara);

deri2 = zeros(ParaNum, ROISize*ROISize);
deri2_ln = zeros(ParaNum, ROISize*ROISize);


%%

for i = 1:ParaNum
    %
    CurInputPara = InputPara;
    fluo = GaussianPSF_w(CurInputPara, ROISize);
    fluo = fluo';
    
    Idat_0dh = fluo(:)';
    Idat_0dh_ln = log(fluo(:))';

    %
    CurInputPara(i)=InputPara(i)+dh;
    fluo = GaussianPSF_w(CurInputPara, ROISize);
    fluo = fluo';

    Idat_pdh = fluo(:)';
    Idat_pdh_ln = log(fluo(:))';

    %
    CurInputPara(i) = InputPara(i)-dh;
    fluo = GaussianPSF_w(CurInputPara, ROISize);
    fluo = fluo';
    Idat_ndh = fluo(:)';
    Idat_ndh_ln = log(fluo(:))';

    % 
    deri2(i,:) = (Idat_pdh - 2*Idat_0dh + Idat_ndh)/dh/dh;
    deri2_ln(i,:) = (Idat_pdh_ln - 2*Idat_0dh_ln + Idat_ndh_ln)/dh/dh;

end



