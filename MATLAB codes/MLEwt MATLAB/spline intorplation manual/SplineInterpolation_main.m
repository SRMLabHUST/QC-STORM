load('testdata.mat')

iData = MeanData2;

DataLen = length(iData);

x = 0:DataLen-1;

y = iData;

plot(x, y,'-o')
hold on

IntGap = 0.1;
xn = 0:IntGap:DataLen+1;

InputDataNum = length(x);

yn = interp1(x, y, xn,'spline');

plot(xn, yn,'-x')

yn1 = SplineInterpolation(x, y, xn);
plot(xn, yn1,'-x')

