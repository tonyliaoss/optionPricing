% clc
% clear all

% set test values...
E = 10;
r = 0.1;
S = 14;
sigma = 0.4;
tau = 0.25;

addpath('simple');
% implicit
coarse = ImplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 1000, ...
	  'numPartitionsT', 1000);
coarse_err = BSEqnEuropean(S, tau, E, r, sigma) - coarse

fine = ImplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 2000, ...
	  'numPartitionsT', 4000);
fine_err = BSEqnEuropean(S, tau, E, r, sigma) - fine
error_ratio = coarse_err/fine_err
fprintf('%f\n%f %f\n',BSEqnEuropean(S, tau, E, r, sigma), coarse, fine);

% fine = ImplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 1000, ...
% 	  'numPartitionsT', 4000);
% fine_err = BSEqnEuropean(S, tau, E, r, sigma) - fine
% error_ratio = coarse_err/fine_err

% Crank-Nicolson
coarse = CrankNicolsonEuropean(S, tau, E, r, sigma, 'numPartitionsX', 1000, ...
	  'numPartitionsT', 1000);
coarse_err = BSEqnEuropean(S, tau, E, r, sigma) - coarse

fine = CrankNicolsonEuropean(S, tau, E, r, sigma, 'numPartitionsX', 2000, ...
	  'numPartitionsT', 2000);
fine_err = BSEqnEuropean(S, tau, E, r, sigma) - fine

error_ratio = coarse_err/fine_err

return

% addpath('sor');
% % implicit
% coarse = sorImplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 1000, ...
% 	  'numPartitionsT', 1000);
% coarse_err = BSEqnEuropean(S, tau, E, r, sigma) - coarse
%
% fine = sorImplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 2000, ...
% 	  'numPartitionsT', 2000);
% fine_err = BSEqnEuropean(S, tau, E, r, sigma) - fine
%
% error_ratio = coarse_err/fine_err
%
% % Crank-Nicolson
% coarse = sorCrankNicolsonEuropean(S, tau, E, r, sigma, 'numPartitionsX', 1000, ...
% 	  'numPartitionsT', 1000);
% coarse_err = BSEqnEuropean(S, tau, E, r, sigma) - coarse
%
% fine = sorCrankNicolsonEuropean(S, tau, E, r, sigma, 'numPartitionsX', 2000, ...
% 	  'numPartitionsT', 2000);
% fine_err = BSEqnEuropean(S, tau, E, r, sigma) - fine
%
% error_ratio = coarse_err/fine_err
%
