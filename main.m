clc
clear all

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

fprintf('  Calculating the closed-form solution and explicit Euler solutions...\n');
% evaluating "simple" numerical methods
bse = BSEqnEuropean(S, tau, E, r, sigma);
addpath('simple');
ee = ExplicitEuropean(S, tau, E, r, sigma);
fprintf('  Using builtin matrix inversion for solving linear systems...\n');
tic
ie = ImplicitEuropean(S, tau, E, r, sigma);
cne = CrankNicolsonEuropean(S, tau, E, r, sigma);
fprintf('    ');
toc

fprintf('  Using SOR for solving linear systems...\n');
% evaluating European options with successive over-relaxation...
addpath('sor');
sor_ee = NaN; % doesn't apply here. SOR only solves linear systems.
tic
sor_ie = SORImplicitEuropean(S, tau, E, r, sigma);
sor_cne = SORCrankNicolsonEuropean(S, tau, E, r, sigma);
fprintf('    ');
toc

% Generate a table.
methodNames = {'Black-Scholes'; 'Explicit'; 'Implicit'; 'Crank-Nicolson'};
simple = [bse; ee; ie; cne];
sor = [bse; sor_ee; sor_ie; sor_cne];
euroTable = table(simple, sor, 'RowNames', methodNames)

