Sel=11;

TotalPn=LocArry(Sel,7);
x0=0;
y0=0;
PSFSigma=LocArry(Sel,5);
Bkg=LocArry(Sel,8);

QE=0.72;
PixelSize=100;

InputPara = [TotalPn*QE x0 y0 PSFSigma Bkg*QE];
ROISize=9;

    %%
    %VarStd is std of each input parameter
    [ParaStd, ~] = GetCRLB(InputPara, ROISize)

    LocPrec=ParaStd(1)*PixelSize

