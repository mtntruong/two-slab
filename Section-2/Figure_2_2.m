% Application of the cut-off projection to solve a backward
% heat conduction problem in a two-slab composite system
% NH Tuan, VA Khoa, MTN Truong, TT Hung, MN Minh

% Copyright (C) 2019

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% Proposed method for finding numerical eigenvalues

clear all; close all; clc;

% Number of the eigenvalues
N = 50;

% a and b as in (−b < x < 0) and (0 < x < a)
a = 3; b = 5;

% For figure 2.2a, K = [K_b, K_a]
K = [1,1]; k_alpha = [1,1];
[ found_roots, lambda_alpha ] = num_eigen( N, a, b, K, k_alpha );

% Plot results:
figure(1)
plot(0:N,found_roots,'ko',0:N,lambda_alpha,'k*')
grid on
xlabel('$N$','Interpreter','latex')
ylabel('$\lambda_{\alpha}$','Interpreter','latex')
h = legend('$\lambda_{b}$','$\lambda_{a}$','Location','northwest');
set(h,'Interpreter','latex','Fontsize',12)

% For figure 2.2b, K = [K_b, K_a]
K = [3.42,1.05]; k_alpha = [0.838,0.339];
[ found_roots, lambda_alpha ] = num_eigen( N, a, b, K, k_alpha );

% Plot results:
figure(2)
plot(0:N,found_roots,'ko',0:N,lambda_alpha,'k*')
grid on
xlabel('$N$','Interpreter','latex')
ylabel('$\lambda_{\alpha}$','Interpreter','latex')
h = legend('$\lambda_{b}$','$\lambda_{a}$','Location','northwest');
set(h,'Interpreter','latex','Fontsize',12)
