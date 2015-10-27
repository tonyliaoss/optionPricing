% set test values...
E = 105;
r = 0.02;
S = 100;
sigma = 0.2;
tau = 0.25;

% print initial conditions to console...
formatString = ['Asset price = %.2f, Strike price = %.2f,' ...
        ' Interest rate = %.4f, Volatility = %.4f, Time to expiry = %.4f\n'];
fprintf(formatString, S, E, r, sigma, tau);

% start off by evaluating European options
fprintf('Evaluating pricing of European options...\n');

% evaluating "simple" numerical methods
bse = BSEqnEuropean(S, tau, E, r, sigma);
addpath('simple');
ee = ExplicitEuropean(S, tau, E, r, sigma);
ie = ImplicitEuropean(S, tau, E, r, sigma);
cne = CrankNicolsonEuropean(S, tau, E, r, sigma);

% evaluating European options with successive over-relaxation...
addpath('sor');
sor_ee = NaN; % doesn't apply here. SOR only solves linear systems.
sor_ie = ImplicitEuropean(S, tau, E, r, sigma);
sor_cne = CrankNicolsonEuropean(S, tau, E, r, sigma);

% Generate a table.
methodNames = {'Black-Scholes'; 'Explicit'; 'Implicit'; 'Crank-Nicolson'};
simple = [bse; ee; ie; cne]
sor = [bse; sor_ee; sor_ie; sor_cne]
tbl = table(simple, sor, 'RowNames', methodNames)

