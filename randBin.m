function [sN] = randBin(N,p, stoFlag)
%randBin draw vectorized random binomial. Length of N and p needs to be
%equal, or one need to be length 1.
%   calls randCppFun.cpp with correct arguments.
if stoFlag == 1
    sN = randCppFun( 0, N, p );
else
    x = reshape(N,prod(size(N)),1);
    y = reshape(p,prod(size(p)),1);
    sN = x .* y;
end
if( (length(N) > 1) || (length(p)==1) )
    sN = reshape(sN, size(N) );
else
    sN = reshape(sN, size(p) );
end
end

