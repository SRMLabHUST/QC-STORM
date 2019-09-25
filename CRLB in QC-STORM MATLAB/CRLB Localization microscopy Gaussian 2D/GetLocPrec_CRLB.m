function CRLB_X = GetLocPrec_CRLB(TotalPhoton, BackgroundPhoton, PSFSigma, PixelSize, QE, ReadNoise, TakeReadNoiseToBackground, IsEMCCD)

% TotalPhoton photon
% BackgroundPhoton photon
% PSFSigma pixel
% PixelSize nm

x0 = 0.0;
y0 = 0.0;
ROISize = 9;


% ReadNoise=1.3; %e-
% TakeReadNoiseToBackground = 1;
% IsEMCCD = 0;

Datlen = length(TotalPhoton);

CRLB_X = zeros(Datlen, 1);

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



%%
for i=1:Datlen
    
    InputPara = [TotalPn_e(i) x0 y0 PSFSigma(i) Bkg_e(i)];

    %%
    %VarStd is std of each input parameter
    [ParaStd, SNR] = GetCRLB(InputPara, ROISize);

    StdX0 = ParaStd(2) * PixelSize;

    CRLB_X(i) = StdX0 ;
end

