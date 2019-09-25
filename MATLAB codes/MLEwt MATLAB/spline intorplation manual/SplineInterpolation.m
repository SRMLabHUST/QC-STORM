function yn = SplineInterpolation(x, y, xn)

InputDataNum = length(x);

%% parameters of spline interpolation

An = y;
Bn = zeros(InputDataNum, 1);
% Cn = zeros(InputDataNum, 1);
Dn = zeros(InputDataNum, 1);


%% construct linear equations to solve Cn

A = zeros(InputDataNum,InputDataNum);
A(1,1) = 1;
A(InputDataNum,InputDataNum) = 1;

for i = 2:InputDataNum-1
    
    h0 = x(i)-x(i-1);
    h1 = x(i+1)-x(i);
    a_t = [h0, 2*(h0+h1), h1];
    A(i,i-1:i+1) = a_t;
    
end

% construct b
b = zeros(InputDataNum, 1);
b(1) = 0;
b(InputDataNum) = 0;

for i = 2:InputDataNum-1
    h0 = x(i) - x(i-1);
    h1 = x(i+1) - x(i);
    
    a0 = y(i-1);
    a1 = y(i);
    a2 = y(i+1);
    
    b(i) = 3/h1*(a2-a1)-3/h0*(a1-a0);
end

%% solve Cn
% Cn = linsolve(A,b)
Cn = SolveSplineEquations(A, b, InputDataNum);

%% solve Bn and Dn
for i = 1:InputDataNum-1
    hi = x(i+1)-x(i);
    Bn(i) = 1/hi*(An(i+1) - An(i)) - hi/3*(2*Cn(i) + Cn(i+1));
    Dn(i) = (Cn(i+1) - Cn(i))/(3*hi);
end

%%

yn = zeros(1,length(xn));

for i = 1:length(xn)
    curx = xn(i);
    
    pos = curx>=x;
    pos = strfind(pos,1);
    if(isempty(pos))
        ParaSet=1;
    else
        ParaSet = pos(end);
    end
    
    
    if(ParaSet<1)
        ParaSet=1;
    end
    if(ParaSet>InputDataNum-1)
        ParaSet=InputDataNum-1;
    end
    
    xi = x(ParaSet);
    
    yn(i) = An(ParaSet) + Bn(ParaSet)*(curx-xi) + Cn(ParaSet)*(curx-xi)^2 + Dn(ParaSet)*(curx-xi)^3;
end



