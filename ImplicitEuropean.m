function  [ P ] = ImplicitEuropean( S, tau, E, r, sigma )
% ImplicitEuropean: Computes the fair value of the European put option
% using the implicit Euler method.
%   Detailed explanation goes here
numPartitionsX = 500;
numPartitionsT = 500;

% determine the maximum and minimum x values.
z = 3; % variable for random walk
S_t = S * exp((r - 0.5 * sigma^2) * tau + sigma * sqrt(tau) * z);
x_max = log(S_t / S);
x_min = -x_max;

% compute the partition size dt and dx.
dt = tau / numPartitionsT;
dx = (x_max - x_min) / numPartitionsX;

% compute coefficients A, B, and C.
A_scalar = (r * dt / (2 * dx) - sigma ^ 2 * dt / (2 * dx ^ 2));
A = ones(numPartitionsX - 1, 1);
A = A * A_scalar;
B_scalar = (sigma ^ 2 * dt / dx ^ 2 + r * dt + 1);
B = ones(numPartitionsX - 1, 1);
B = B * B_scalar;
C_scalar = (-(r * dt) / (2 * dx) - (sigma ^ 2 * dt) / (2 * dx ^ 2));
C = ones(numPartitionsX - 1, 1);
C = C * C_scalar;

% generate a tri-diagonal matrix T, of size numPartitionsX - 1.
TRI = spdiags([A, B, C], [-1, 0, 1], numPartitionsX - 1, numPartitionsX - 1);

% allocate memory for the solution mesh
PRICE = zeros(numPartitionsT + 1, numPartitionsX + 1);
X = linspace(x_min, x_max, numPartitionsX + 1);
T = linspace(0, tau, numPartitionsT + 1);

% initial condition...
PRICE(1, :) = max(E - S * exp(X), 0);

% boundary value at X = -inf
PRICE(:, 1) = E * exp(-r * T); % present value of payoff
% boundary value at X = inf
PRICE(:, end) = 0;

P_boundary = zeros(numPartitionsX - 1, 1);
% populate the solution mesh
for i = 2:numPartitionsT + 1
    P_boundary(1) = A_scalar * PRICE(i, 1);
    P_boundary(end) = C_scalar * PRICE(i, end);
    % we want to solve the equation: TRI * P(t = t) = P(t = t - 1) + P_boundary(t = t)
    PRICE(i, 2:numPartitionsX) = ...
        transpose((transpose(PRICE(i-1, 2:numPartitionsX)) - P_boundary) \ TRI);
end

P = PRICE(numPartitionsT+1, ceil((numPartitionsX + 1) / 2));

end

