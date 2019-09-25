function ParaStd_all = GetLocPrec_CRLB3D(TotalPhoton, BackgroundPhoton, PSFSigmaX, PSFSigmaY, QE, ReadNoise, TakeReadNoiseToBackground, IsEMCCD, ROISize)

x0 = 0.0;
y0 = 0.0;


TotalFluoNum = length(TotalPhoton);

CRLB_XY_nm = zeros(TotalFluoNum, 2);

% shot noise=sqrt(pn*QE)

%%
if(IsEMCCD==0)
    % CMOS camera
    QE_e = QE;
    if(TakeReadNoiseToBackground == 1)
        BackgroundPhoton = (BackgroundPhoton*QE + ReadNoise*ReadNoise)/QE;
    end

else
    QE_e = QE/sqrt(2);
    if(TakeReadNoiseToBackground == 1)
        BackgroundPhoton = (BackgroundPhoton*2*QE + ReadNoise*ReadNoise)/(2*QE);
    end
end

TotalPn_e = TotalPhoton*QE_e;
Bkg_e = BackgroundPhoton*QE_e;

ParaStd_all = zeros(TotalFluoNum,6);

for i=1:TotalFluoNum
    if(mod(i,round(TotalFluoNum/20))==0)
        disp(i)
    end
    InputPara = [TotalPn_e(i) x0 y0 PSFSigmaX(i) PSFSigmaY(i) Bkg_e(i)];
    
    ParaStd = GetCRLB(InputPara, ROISize);
    
    ParaStd_all(i,:) = ParaStd;
end






