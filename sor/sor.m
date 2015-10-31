function  [ x ] = sor( A, b, x0 )
% sor: Solves a system of linear equations using the successive over-relaxation
% technique. Note that the error tolerance and maximum number of iterations are
% hard-coded into this function.
% Solves the linear system Ax = b.
%
%   Detailed explanation goes here
% INPUT PARAMETERS
%   A - the input matrix corresponding to a linear transformation in the
%       equation Ax = b, where x is the unknown.
%   b - the input column vector corresponding to the right hand side of the
%       equation Ax = b, where x is the unknown.
%   x0 - the initial "guess" at the vector x.
% OUTPUT PARAMETERS
%   x - the output column vector corresponding to the unknown vector in the
%       equation Ax = b.

% hard-coded parameters.
% relaxation rate. Method becomes Gauss-Seidel method when omega = 1.
omega = 1.5; % empirically, seems to be the best number.
% % maximum number of iterations
% maxIter = 2000;
% error tolerance
epsilon = 1e-8;

% split the matrix into 3 parts...
[L, D, U] = ldu(A);

% create a local copy of the inital guess
x = x0;
% make xold big so that we enter the loop below
xold = 16384 * x;

numIter = 0;
while norm(xold - x) > epsilon
  % this corresponds to the 'previous' guess
  xold = x;

%  % iterate through each component of x... We need to do this with a for loop
%  % because we want to use updated values as soon as it becomes available.
%  % i.e. use Gauss-Seidel formulation instead of Jacobi formulation.
%  % matlab doesn't like for loops... but in real, fast code written in C++
%  this is the approach to try first.
%  for i = 1:length(x)
%    x(i) = (1 - omega) * x(i)...
%	    + omega / D(i,i) * (b(i) - L(i,:) * x - U(i,:) * x);
%  end

  % now we update x according to the formula found on Wikipedia :p
  % This is the Jacobi iteration... where we compute *all* of x for the
  % iteration k before proceeding to iteration k+1. Theoretically the
  % Gauss-Seidel formulation should be faster... but I guess Matlab is just
  % insanely fast with arrays...? :s
  x = (D + omega * L) \ ...
        (omega * b - (omega * U + (omega - 1) * D) * x);

  numIter = numIter + 1;
end

end

function [ L, D, U ] = ldu( A )
% ldu: decomposes the square matrix A into its strictly-lower-diagonal
% component L, its strictly-upper-diagonal component U, and a diagonal matrix
% D. Note that this is different from LU factorization. LU factorization
% requires that L * U = A. In here we need L + D + U = A instead.
%   Detailed explanation goes here
% INPUT PARAMETERS
%   A - the input square matrix
% OUTPUT PARAMETERS
%   L - the strictly lower triangular matrix with 0 diagonal terms.
%   D - the strictly diagonal matrix with nonzero diagonal terms and zeros
%       everywhere else.
%   U - the strictly upper triangular matrix with 0 diagonal terms.

L = tril(A, -1);
D = diag(diag(A));
U = triu(A, 1);

end

