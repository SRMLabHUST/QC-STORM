function ZCrlb = GetZLocPrec_all(PSFSigmaX, PSFSigmaY, LocPrec_X, LocPrec_Y, PixelSize, Calib3d_fitresult)

FluoNum = length(PSFSigmaX);
ZCrlb = zeros(FluoNum,1);

for i = 1:FluoNum
    if(mod(i,round(FluoNum/20))==0)
       disp(i); 
    end
    
    PSFSigmaX_t = PSFSigmaX(i);
    PSFSigmaY_t = PSFSigmaY(i);
    LocPrec_X_t = LocPrec_X(i);
    LocPrec_Y_t = LocPrec_Y(i);
    ZCrlb(i) = GetZLocPrec(PSFSigmaX_t, PSFSigmaY_t, LocPrec_X_t, LocPrec_Y_t, PixelSize, Calib3d_fitresult);
end

