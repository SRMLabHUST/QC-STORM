function [X0,Y0,SigmaX,SigmaY] = SigmaPreEstimate(moleculeSub)

ImgDat=double(moleculeSub);
%%
mx1=mean(ImgDat,1);
my1=mean(ImgDat,2);

% mx1=mx1-(mx1(1)+mx1(end))/2;
% my1=my1-(my1(1)+my1(end))/2;

% plot(mx1)
% hold on
% plot(my1)

DatLen=length(mx1);

[MaxDatX,MaxPosX]=max(mx1);
[MaxDatY,MaxPosY]=max(my1);

X0=MaxPosX-0.5;
Y0=MaxPosY-0.5;
%%
WDat=0;
SDat=0;

for i=1:MaxPosX-1
    WDat = WDat + mx1(i)*abs(i-MaxPosX);
    SDat = SDat + mx1(i);
    
end
SigmaL=WDat/SDat;
WDat=0;
SDat=0;

for i = MaxPosX+1:DatLen
    WDat = WDat + mx1(i)*abs(i-MaxPosX);
    SDat = SDat + mx1(i);
    
end
SigmaR=WDat/SDat;

SigmaX=(SigmaL+SigmaR)/2;

%%
WDat=0;
SDat=0;

for i=1:MaxPosY-1
    WDat = WDat + my1(i)*abs(i-MaxPosY);
    SDat = SDat + my1(i);
    
end
SigmaA=WDat/SDat;
WDat=0;
SDat=0;

for i = MaxPosY+1:DatLen
    WDat = WDat + my1(i)*abs(i-MaxPosY);
    SDat = SDat + my1(i);
    
end
SigmaB=WDat/SDat;

SigmaY=(SigmaA+SigmaB)/2;
