function [ P ] = montecarloAmerican( S, tau, E, r, sigma , varargin )
% Gives the price of an American put option without dividends, as computed
% through Monte-Carlo simulation by least squares (Longstaff and Schwartz).
%   Detailed explanation goes here
% INPUT PARAMETERS
%   S - the current market value of the asset.
%   tau - time to expiry of the option, expressed in years.
%   E - the strike price, or exercise price, of the option.
%   r - the risk-free investment return per annum.
%   sigma - the market volatility.

%% note also need to choose number and type of basis functions.
% e.g. number = 2 and type = polynomial -> least squares.

% default values for optional arguments...
defaultN = 500; % number of sample paths.
defaultNumPartitionsT = 600; % number of observation dates.

% define inputParser to parse optional arguments.
parser = inputParser;
parser.addRequired('S', @isnumeric);
parser.addRequired('tau', @isnumeric);
parser.addRequired('E', @isnumeric);
parser.addRequired('r', @isnumeric);
parser.addRequired('sigma', @isnumeric);
parser.addParamValue('N', defaultN, @isnumeric);
parser.addParamValue('numPartitionsT', defaultNumPartitionsT, @isnumeric);

% parse input and setting optional values...
parser.parse(S, tau, E, r, sigma, varargin{:});
inputs = parser.Results;
N = inputs.N;
numPartitionsT = inputs.numPartitionsT;

dt = tau / numPartitionsT;

% generate N paths and store them all.
sim = simulateBrownian(S, r, sigma, dt, numPartitionsT, N);

% and then run the regression.
sim(:, end) = max(0, sim(:, 1) - sim(:, end));
sim
A = [];
for i = numPartitionsT:-1:1
  % consider only in-the-money options (i) < (1)
  % W is the cash flow realized if we exercise right now.
  W = max(0, sim(:, 1) - sim(:, i));
  % X is the current price of asset.
  X = sim(:, i);
  % Y is the discounted cash flow if we exercise later.
  Y = exp(-r * dt) * sim(:, i+1);
%  Y(X==0) = 0; % ignore any out-of-the-money options

  % generate two "reduced" vectors so that we can do regression
  X_red = X;
  Y_red = Y;
  X_red(W==0) = []; % ignore any out-of-money options
  Y_red(W==0) = [];

  % regress Y on X.
  coeff = lseRegression(X_red, Y_red);
  % unpack the coefficients...
  a = coeff(1);
  b = coeff(2);
  c = coeff(3);

  % Z is the expected cash flow from continuing the option conditional on the
  % time step right now.
  regressed_fn = @(X, W) (W ~= 0) * (a * X .^ 2 + b * X + c); % returns 0 if X == 0
  % calculate the conditional expectation according to regressed parameters
  Z = arrayfun(regressed_fn, X, W);

  % the next step is to compare W and Z I guess?
  % if W <= Z, keep. (do nothing)
  % if W > Z, exercise.
  exercise = zeros(N, 1); % a "vectorized" boolean variable.
  exercise(W > Z) = 1;
  % if don't exercise, discount future cash flow, else update the value to W.
  sim(:, i) = (1 - exercise) .* Y + exercise .* W;
  % if we exercise now, then we can't exercise later by definition
  sim(:, i+1) = zeros(N, 1);
%  sim(:, i+1) = (1 - exercise) .* sim(:, i+1) + exercise .* 0;
end
sim

% count only the exercised paths.
P = mean(sim(:, 1));

end

