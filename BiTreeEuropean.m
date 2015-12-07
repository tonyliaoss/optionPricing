function [ output_args ] = BiTreeEuropean( S, tau, E, r, sigma, varargin )
% BiTreeEuropean: Computes the fair value of the European put option
% using the binomial tree method
% INPUT PARAMETERS
%   S - the current market value of the asset.
%   tau - time to expiry of the option, expressed in years.
%   E - the strike price, or exercise price, of the option.
%   r - the risk-free investment return per annum.
%   sigma - the market volatility.
%   varargin - optional arguments as explained below.
% OPTIONAL INPUT (PARAMETER-VALUE PAIRS)
%   'numPartitionsT', nT - number of partitions in T

% default values for optional arguments...
defaultNumPartitionsT = 128; % days in a year

% define inputParser to parse optional arguments.
parser = inputParser;
parser.addRequired('S', @isnumeric);
parser.addRequired('tau', @isnumeric);
parser.addRequired('E', @isnumeric);
parser.addRequired('r', @isnumeric);
parser.addRequired('sigma', @isnumeric);
parser.addParamValue('numPartitionsT', defaultNumPartitionsT, @isnumeric);

% parse input and setting optional values...
parser.parse(S, tau, E, r, sigma, varargin{:});
inputs = parser.Results;
numPartitionsT = inputs.numPartitionsT;

% binomial parameters... (from Prentice Hall)
u = exp(sigma * sqrt(tau)); % proprotion of increase for moving up
d = 1 / u; % proportion of decrease for moving down

dt = tau / numPartitionsT;
a = exp(r * sqrt(dt));
p = (a - d) / (u - d); % "probability" of up move

% todo

end

