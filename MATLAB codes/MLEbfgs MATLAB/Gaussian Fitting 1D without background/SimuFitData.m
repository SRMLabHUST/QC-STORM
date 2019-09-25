

Amp0=25806;
X0=51;
Sigma0= 3.8;

Amp1=4886;
X1=61.8;
Sigma1=12;

DatLen=150;
Sigma01=1/(2*Sigma0*Sigma0);
Sigma11=1/(2*Sigma1*Sigma1);

XPos = (1:1:DatLen) - 0.5;
YDat0 = Amp0 * exp(-(XPos - X0).^2 * Sigma01);
YDat1 = Amp1 * exp(-(XPos - X1).^2 * Sigma11);
YDat=int32(YDat0+YDat1);

plot(XPos,YDat)


save stdimg XPos YDat Amp0 X0 Sigma0

%%

WriteDat=int32(y);

FileName='fitArry_whole.txt';

WriteDat=WriteDat';
WriteDat=int32(WriteDat);

fid=fopen(FileName,'w');
fwrite(fid,WriteDat(:),'int32');

fclose(fid);


