clc
clear all

% set test values...
E = 105;
r = 0.02;
S = 100;
sigma = 0.17;
tau = 0.25;

% sweep through asset prices to get a plot of convergence rates...
addpath('simple');
% % explicit
% asset_prices = [];
% error_ratios = [];
% figure
% for S=8:0.5:12
%   coarse = ExplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 200);
%   coarse_err = BSEqnEuropean(S, tau, E, r, sigma) - coarse;
%   fine = ExplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 400);
%   fine_err = BSEqnEuropean(S, tau, E, r, sigma) - fine;
%   error_ratio = coarse_err/fine_err;
%
%   asset_prices = [asset_prices S];
%   error_ratios = [error_ratios error_ratio];
% end
% scatter(asset_prices, error_ratios);
% a = error_ratios';
% b = num2str(a);
% c = cellstr(b);
% dx = 0.1; dy = 0.1;
% text(asset_prices + dx, error_ratios + dy, c);
% title('Error ratios as a function of initial asset prices (Explicit European).');
% xlabel('Initial asset price');
% ylabel('Error ratio');
%
% % implicit
% asset_prices = [];
% error_ratios = [];
% figure
% for S=8:0.5:12
%   coarse = ImplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 1000, ...
%   	  'numPartitionsT', 1000);
%   coarse_err = BSEqnEuropean(S, tau, E, r, sigma) - coarse;
%   fine = ImplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 2000, ...
%   	  'numPartitionsT', 4000);
%   fine_err = BSEqnEuropean(S, tau, E, r, sigma) - fine;
%   error_ratio = coarse_err/fine_err;
%
%   asset_prices = [asset_prices S];
%   error_ratios = [error_ratios error_ratio];
% end
% scatter(asset_prices, error_ratios);
% a = error_ratios';
% b = num2str(a);
% c = cellstr(b);
% dx = 0.1; dy = 0.1;
% text(asset_prices + dx, error_ratios + dy, c);
% title('Error ratios as a function of initial asset prices (Implicit European).');
% xlabel('Initial asset price');
% ylabel('Error ratio');

% crank-nicolson
asset_prices = [];
error_ratios = [];
hold on
for S=100:0.2:110
  coarse = CrankNicolsonEuropean(S, tau, E, r, sigma, 'numPartitionsX', 500, ...
  	  'numPartitionsT', 500);
  coarse_err = BSEqnEuropean(S, tau, E, r, sigma) - coarse;
  fine = CrankNicolsonEuropean(S, tau, E, r, sigma, 'numPartitionsX', 1000, ...
  	  'numPartitionsT', 1000);
  fine_err = BSEqnEuropean(S, tau, E, r, sigma) - fine;
  error_ratio = coarse_err/fine_err;

  asset_prices = [asset_prices S];
  error_ratios = [error_ratios error_ratio];
end
% plot(asset_prices, error_ratios);
scatter(asset_prices, error_ratios);
a = error_ratios';
b = num2str(a);
c = cellstr(b);
dx = 0.05; dy = 0.05;
text(asset_prices + dx, error_ratios + dy, c);
title('Error ratios as a function of initial asset prices (Crank-Nicolson European).');
xlabel('Initial asset price');
ylabel('Error ratio');

% set test values...
E = 10;
r = 0.1;
S = 14;
sigma = 0.4;
tau = 0.25;

return % stop here

%% explicit
%strike_prices = [];
%error_ratios = [];
%figure
%for E=0:2:14
%  coarse = ExplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 200);
%  coarse_err = BSEqnEuropean(S, tau, E, r, sigma) - coarse;
%  fine = ExplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 400);
%  fine_err = BSEqnEuropean(S, tau, E, r, sigma) - fine;
%  error_ratio = coarse_err/fine_err;
%
%  strike_prices = [strike_prices E];
%  error_ratios = [error_ratios error_ratio];
%end
%scatter(strike_prices, error_ratios);
%a = error_ratios';
%b = num2str(a);
%c = cellstr(b);
%dx = 0.1; dy = 0.1;
%text(asset_prices + dx, error_ratios + dy, c);
%title('Error ratios as a function of strike prices (Explicit European).');
%xlabel('Strike price');
%ylabel('Error ratio');
%
%% Implicit
%strike_prices = [];
%error_ratios = [];
%figure
%for E=0:2:14
%  coarse = ImplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 400, ...
%  	  'numPartitionsT', 400);
%  coarse_err = BSEqnEuropean(S, tau, E, r, sigma) - coarse;
%  fine = ImplicitEuropean(S, tau, E, r, sigma, 'numPartitionsX', 800, ...
%  	  'numPartitionsT', 800);
%  fine_err = BSEqnEuropean(S, tau, E, r, sigma) - fine;
%  error_ratio = coarse_err/fine_err;
%
%  strike_prices = [strike_prices E];
%  error_ratios = [error_ratios error_ratio];
%end
%scatter(strike_prices, error_ratios);
%a = error_ratios';
%b = num2str(a);
%c = cellstr(b);
%dx = 0.1; dy = 0.1;
%text(asset_prices + dx, error_ratios + dy, c);
%title('Error ratios as a function of strike prices (Implicit European).');
%xlabel('Strike price');
%ylabel('Error ratio');
%
%% crank-nicolson
%strike_prices = [];
%error_ratios = [];
%figure
%for E=0:2:14
%  coarse = CrankNicolsonEuropean(S, tau, E, r, sigma, 'numPartitionsX', 1000, ...
%  	  'numPartitionsT', 1000);
%  coarse_err = BSEqnEuropean(S, tau, E, r, sigma) - coarse;
%  fine = CrankNicolsonEuropean(S, tau, E, r, sigma, 'numPartitionsX', 2000, ...
%  	  'numPartitionsT', 2000);
%  fine_err = BSEqnEuropean(S, tau, E, r, sigma) - fine;
%  error_ratio = coarse_err/fine_err;
%
%  strike_prices = [strike_prices E];
%  error_ratios = [error_ratios error_ratio];
%end
%scatter(strike_prices, error_ratios);
%a = error_ratios';
%b = num2str(a);
%c = cellstr(b);
%dx = 0.1; dy = 0.1;
%text(asset_prices + dx, error_ratios + dy, c);
%title('Error ratios as a function of strike prices (Crank-Nicolson European).');
%xlabel('Strike price');
%ylabel('Error ratio');
%
