function [ P ] = BiTreeEuropean( S, tau, E, r, sigma, varargin )
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
%   'numPartitionsT', nT - number of time steps.

% default values for optional arguments...
defaultNumPartitionsT = 256;

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

%%%% binomial parameters...
dt = tau / numPartitionsT;

%% alternate way to set up the parameters.
% A = 0.5 * (exp(-r* dt) + exp((r + sigma ^ 2) * dt));
%
% u = A + sqrt(A^2 - 1);
% d = 1/u;

u = exp(sigma * sqrt(dt));
d = 1/u;
a = exp(r * dt);
p = (a - d) / (u - d); % "probability" of up move

% how to build a tree? We initialize a vector of length 1, and gradually add to
% it...
Tree = S; % initialization

for i = 1:numPartitionsT
  mult = [u*ones(1,i) d]; % [u u ... u d]
  Tree = [Tree Tree(i)];

  Tree = Tree .* mult;
end
Tree = Tree'; % make sure Tree is column

Tree = max(E - Tree, 0); % option value at expiration date.

% backpropagation matrix
B = spdiags([p*ones(numPartitionsT+1, 1), (1-p)*ones(numPartitionsT+1, 1)], ...
    [0, 1], zeros(numPartitionsT+1, numPartitionsT+1));

for i = numPartitionsT:-1:1
  B(1:i, 1:i+1);
  Tree = B(1:i, 1:i+1) * Tree;
end

P = exp(-r * tau) * Tree;

end

