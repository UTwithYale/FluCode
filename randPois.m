function [sN] = randPois(mu, stoFlag)
%randBin draw vectorized random poissson. 
%   calls randCppFun.cpp with correct arguments.
if stoFlag == 1
    sN = reshape( randCppFun( 2, mu, 0 ), size(mu) );
else
    sN = mu;
end
end

