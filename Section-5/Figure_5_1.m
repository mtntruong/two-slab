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

clear all; close all; clc;

% Parameters

epsilon = 10^-6;

a = 1; b = 1; c = 1;

kappaA = 0.339;
kappaB = 0.838;

kA = 1.05;
kB = 3.42;

t0 = 0; tf = 0.1;

M = 30; N = 30;

% Initialize

load lambdas.mat
% You can re-produce lambdas.mat by running the following script
%[lambda_bar, nu_a, nu_b] = newThetaEps( M, N, a, b, c, t0, tf, kA, kB, kappaA, kappaB, epsilon );


[X,Y] = meshgrid(0:30,0:30);
l_a = []; l_b = [];
for m = 0 : M
    for n = 0 : N
        l_a(m+1,n+1) = nu_a(m+1,n+1)^2 + ((m)*pi)^2;
        l_b(m+1,n+1) = nu_b(m+1,n+1)^2 + ((m)*pi)^2;
    end
end

figure(1)
hSurf = surf(X,Y,l_a);
xlabel('$n$','Interpreter','latex'); ylabel('$m$','Interpreter','latex'); zlabel('$\lambda_{a}^{2}$','Interpreter','latex');
colormap gray

figure(2)
gSurf = surf(X,Y,l_b);
xlabel('$n$','Interpreter','latex'); ylabel('$m$','Interpreter','latex'); zlabel('$\lambda_{b}^{2}$','Interpreter','latex');
colormap gray
