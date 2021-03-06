function  [ P ] = ImplicitEuropean( S, tau, E, r, sigma, varargin )
% ImplicitEuropean: Computes the fair value of the European put option
% using the implicit Euler method.
% INPUT PARAMETERS
%   S - the current market value of the asset.
%   tau - time to expiry of the option, expressed in years.
%   E - the strike price, or exercise price, of the option.
%   r - the risk-free investment return per annum.
%   sigma - the market volatility.
%   varargin - optional arguments as explained below.
% OPTIONAL INPUT (PARAMETER-VALUE PAIRS)
%   'numPartitionsX', nX - number of partitions in X
%   'numPartitionsT', nT - number of partitions in T

% default values for optional arguments...
defaultNumPartitionsX = 2000;
defaultNumPartitionsT = 2000;

% define inputParser to parse optional arguments.
parser = inputParser;
parser.addRequired('S', @isnumeric);
parser.addRequired('tau', @isnumeric);
parser.addRequired('E', @isnumeric);
parser.addRequired('r', @isnumeric);
parser.addRequired('sigma', @isnumeric);
parser.addParamValue('numPartitionsX', defaultNumPartitionsX, @isnumeric);
parser.addParamValue('numPartitionsT', defaultNumPartitionsT, @isnumeric);

% parse input and setting optional values...
parser.parse(S, tau, E, r, sigma, varargin{:});
inputs = parser.Results;
numPartitionsX = inputs.numPartitionsX;
numPartitionsT = inputs.numPartitionsT;

% determine the maximum and minimum x values.
z = 5; % variable for random walk
S_t = S * exp((r - 0.5 * sigma^2) * tau + sigma * sqrt(tau) * z);
if S == 0
  x_max = 2 * E;
  x_min = -2 * E;
else
  x_max = log(S_t / S);
  x_min = -x_max;
end

% compute the partition size dt and dx.
dt = tau / numPartitionsT;
dx = (x_max - x_min) / numPartitionsX;

% compute coefficients A, B, and C.
A_scalar = (((r - 0.5 * sigma ^ 2) * dt) / (2 * dx) - sigma ^ 2 * dt / (2 * dx ^ 2));
A = ones(numPartitionsX - 1, 1);
A = A * A_scalar;
B_scalar = (sigma ^ 2 * dt / (dx ^ 2) + r * dt + 1);
B = ones(numPartitionsX - 1, 1);
B = B * B_scalar;
C_scalar = (-((r - 0.5 * sigma ^ 2) * dt) / (2 * dx) - (sigma ^ 2 * dt) / (2 * dx ^ 2));
C = ones(numPartitionsX - 1, 1);
C = C * C_scalar;

% generate a tri-diagonal matrix TRI, of size numPartitionsX - 1.
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

P_boundary = zeros(1, numPartitionsX - 1);
% populate the solution mesh
for i = 2:numPartitionsT + 1
    P_boundary(1) = A_scalar * PRICE(i, 1);
    P_boundary(end) = C_scalar * PRICE(i, end);

    % we want to solve the equation: TRI * P(t = t) = P(t = t - 1) + P_boundary(t = t)
    PRICE(i, 2:numPartitionsX) = ...
        transpose(TRI \ ...
                    transpose(PRICE(i-1, 2:numPartitionsX) ...
                              - P_boundary) ...
                 );
end

% interpolate the price at x = 0...
P = interp1(X, PRICE(end, :), 0);

end

