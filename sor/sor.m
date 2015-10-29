function  [ x ] = sor( A, b )
% sor: Solves a system of linear equations using the successive over-relaxation
% technique. Note that the error tolerance and maximum number of iterations are
% hard-coded into this function.
% Solves the linear system Ax = b.
%
%   Detailed explanation goes here
% INPUT PARAMETERS
%   A = the input matrix corresponding to a linear transformation in the
%       equation Ax = b, where x is the unknown.
%   b = the input column vector corresponding to the right hand side of the
%       equation Ax = b, where x is the unknown.
% OUTPUT PARAMETERS
%   x = the output column vector corresponding to the unknown vector in the
%       equation Ax = b.

% hard-coded parameters.
% relaxation rate
% omega =
% maximum number of iterations
% maxIter =
% error tolerance
% epsilon =

x = A \ b;

end

