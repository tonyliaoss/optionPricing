function  [ P ] = sorCrankNicolsonEuropean( S, tau, E, r, sigma )
% sorCrankNicolsonEuropean: Computes the fair value of the European put option
% using the Crank-Nicolson implicit finite differences scheme. Uses SOR to
% solve the system of linear equations at each time step.
%   Detailed explanation goes here
% INPUT PARAMETERS
%   S - the current market value of the asset.
%   tau - time to expiry of the option, expressed in years.
%   E - the strike price, or exercise price, of the option.
%   r - the risk-free investment return per annum.
%   sigma - the market volatility.
numPartitionsX = 2000;
numPartitionsT = 2000;

% determine the maximum and minimum x values.
z = 5; % variable for random walk
S_t = S * exp((r - 0.5 * sigma^2) * tau + sigma * sqrt(tau) * z);
x_max = log(S_t / S);
x_min = -x_max;

% compute the partition size dt and dx.
dt = tau / numPartitionsT;
dx = (x_max - x_min) / numPartitionsX;

% compute coefficients A, B, C, D, E, and F.
A_scalar = ((r - 0.5 * sigma ^ 2) * dt / (4 * dx) - sigma ^ 2 * dt / (4 * dx ^ 2));
A = ones(numPartitionsX - 1, 1);
A = A * A_scalar;
B_scalar = (0.5 * sigma ^ 2 * dt / dx ^ 2 + 0.5 * r * dt + 1);
B = ones(numPartitionsX - 1, 1);
B = B * B_scalar;
C_scalar = (-((r - 0.5 * sigma ^ 2) * dt) / (4 * dx) - (sigma ^ 2 * dt) / (4 * dx ^ 2));
C = ones(numPartitionsX - 1, 1);
C = C * C_scalar;
alpha_scalar = - A_scalar;
alpha = ones(numPartitionsX - 1, 1);
alpha = alpha * alpha_scalar;
beta_scalar = ( -0.5 * sigma ^ 2 * dt / dx ^ 2 - 0.5 * r * dt + 1);
beta = ones(numPartitionsX - 1, 1);
beta = beta * beta_scalar;
gamma_scalar = - C_scalar;
gamma = ones(numPartitionsX - 1, 1);
gamma = gamma * gamma_scalar;

% for the Crank-Nicolson method we need to generate two tri-diagonal matrices.
% generate a tri-diagonal matrix TRI_t, of size numPartitionsX - 1, at time t.
TRI_t = spdiags([A, B, C], [-1, 0, 1], numPartitionsX - 1, numPartitionsX - 1);
% generate a tri-diagonal matrix TRI_tplus, of size numPartitionsX - 1, at time t - 1.
TRI_tprev = spdiags([alpha, beta, gamma], [-1, 0, 1], numPartitionsX - 1, numPartitionsX - 1);

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

P_boundary = zeros(1, numPartitionsX - 1);
% populate the solution mesh
for i = 2:numPartitionsT + 1
    P_boundary(1) = A_scalar * PRICE(i, 1);
    P_boundary(end) = C_scalar * PRICE(i, end);

    PRICE_prev = transpose(PRICE(i-1, 2:numPartitionsX));
    % this is the 'right hand side' of the system of linear equation that we
    % wanna solve.
    rhs = (TRI_tprev * PRICE_prev ...
             - transpose(P_boundary));
    % we want to solve the equation:
    % TRI_t * P(t = t) = TRI_tprev * P(t = t - 1) - P_boundary(t = t)
    PRICE(i, 2:numPartitionsX) = ...
        transpose(sor(TRI_t, rhs, PRICE_prev));
end
P = PRICE(numPartitionsT+1, ceil((numPartitionsX + 1) / 2));

% interpolate the price at x = 0...
P = interp1(X, PRICE(end, :), 0);

end

