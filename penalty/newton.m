function [ P ] = newton(rhs, TRI, P0, G, zeta_func, zeta)
% newton: Solves a system of nonlinear equations using Newton iterations. In
% essence it follows the procedure:
%	x(k+1) = x(k) - F'(x(k)) \ F(x(k))
% where x(k) is the kth iteration of x, and F'(x) is the Jacobian of F.

epsilon = 1e-8;

P = P0;

old_P = P * 10000;

%%G
%size(G)
%%P
%size(P)

% zeta_func = @(G, P) (G > P) * zeta;
[R C] = size(TRI); % get the row and column count for TRI.

while norm(old_P - P) > epsilon
  old_P = P;
  zeta_vec = arrayfun(zeta_func, G, P);

%  P = (TRI + diag(zeta_vec)) \ (P0 + zeta_vec .* G);

  zeta_diag = spdiags(zeta_vec, [0], R, C);
  F = (TRI + zeta_diag) * P - (rhs + zeta_vec .* G);
  F_jac = TRI + zeta_diag; % - diag(G .* zeta_vec); % diff(zeta) = -zeta

  P = P - F_jac \ F;
%  norm(old_P - P)
%  norm(P)
end

% while norm(old_x - x) > epsilon
%   old_x = x;
%   x = x - A_jac \ A;
% end

