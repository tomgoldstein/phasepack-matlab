function [ output_args ] = checkAdjointVerbose( A,At,n)

if isnumeric(A)
    x = randn(n,1);
    Ax = A*x;
    y = randn(numel(Ax),1);
    Aty = At*y;
    
    inner1 = Ax'*y;
    inner2 = x'*Aty;
    
    fprintf('inner1 = %f, inner2 = %f, error = %f\n',inner1,inner2,abs(inner1-inner2)/abs(inner1));
else
    x = randn(n,1);
    Ax = A(x);
    y = randn(numel(Ax),1);
    Aty = At(y);
    
    inner1 = Ax'*y;
    inner2 = x'*Aty;
    
    fprintf('inner1 = %f, inner2 = %f, error = %f\n',inner1,inner2,abs(inner1-inner2)/abs(inner1));
end

end

