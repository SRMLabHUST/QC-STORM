function x = SolveSplineEquations(A, b, InputDataNum)

% Xs = linsolve(A, b);

% the matrix of spline interpolation is special
% the matrix is belt like and many zeros
% LU decomposition

% Ax = b
% LUx = b
% Ly = b
% Ux = y


y = zeros(InputDataNum, 1);
x = zeros(InputDataNum, 1);
Beta = zeros(InputDataNum, 1);

for i = 1:InputDataNum-1
    ci = A(i, i+1);
    bi = A(i, i);

    if(i==1)
        Beta(i) = ci/bi;
    else
        ai = A(i, i-1);
        Beta(i) = ci/(bi - ai*Beta(i-1));
    end
end

for i = 1:InputDataNum
    bi = A(i, i);
    
    if(i==1)
        y(i) = b(i)/bi;
    else
        ai = A(i, i-1);
        y(i) = (b(i) - ai*y(i-1))/(bi-ai*Beta(i-1));
    end
end

for i = InputDataNum:-1:1
    if(i==InputDataNum)
        x(i) = y(i);
    else
        x(i) = y(i) - Beta(i)*x(i+1);
    end
end


