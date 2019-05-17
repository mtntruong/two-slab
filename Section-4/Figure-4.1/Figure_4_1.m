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
epsilon = 10^-2;
a = 3; b = 5; C = 1;
t0 = 0; tf = 0.1;
kappaA = 0.339;  kappaB = 0.838;
kA = 1.05;       kB = 3.42;

% Initialize
[nSet, lambda_bar, lambda_a, lambda_b] = ThetaEps( a, b, t0, tf, kA, kB, kappaA, kappaB, epsilon );
J = 50;
generateCN( J, a, b, C, epsilon, t0, tf, lambda_a, lambda_b );

PTa = sym(0);
PTb = sym(0);

% Symbolic

syms x y t

Tb = C * cos(lambda_b(2) * (x + b)) * exp(-t) / ...
    ((a + b) * cos (lambda_b(2) * b));
Ta = C * cos(lambda_a(2) * (x - a)) * exp(-t) / ...
    ((a + b) * cos (lambda_a(2) * a));

% Iteration

for i = 1 : length(nSet)
    n = nSet(i);

    phiA = cos(lambda_a(n) * (x - a)) / cos(lambda_a(n) * a);
    phiB = cos(lambda_b(n) * (x + b)) / cos(lambda_b(n) * b);

    %-- Ta = subs(y + a, tf)
    u = subs(Ta, {x, t}, {y + a, tf});
    if (lambda_a(n) == 0)
        tmpGrid = linspace(-a, 0, 50);
        filonA = trapz(eval(subs(u + epsilon * noise, tmpGrid)));
    else
        filonA = (1 / cos(lambda_a(n) * a)) * adaptive_filon(y, -a, 0, u, lambda_a(n), epsilon);
    end

    %-- Tb = subs(y - b, tf)
    u = subs(Tb, {x, t}, {y - b, tf});
    if (lambda_b(n) == 0)
        tmpGrid = linspace(0, b, 50);
        filonB = trapz(eval(subs(u + epsilon * noise, tmpGrid)));
    else
        filonB = (1 / cos(lambda_b(n) * b)) * adaptive_filon(y, 0, b, u, lambda_b(n), epsilon);
    end

    C_n_A = C_n_new(n, 'a');
    PTa = PTa + ...
        C_n_A * ...
        exp(lambda_bar(n) * (tf - t)) * ...
        phiA;

    C_n_B = C_n_new(n, 'b');
    PTb = PTb + ...
        C_n_B * ...
        exp(lambda_bar(n) * (tf - t)) * ...
        phiB;
end



% 3D plot

figure(1)
ezsurf(Tb, [0 tf -b 0])
hold on
ezsurf(Ta, [0 tf 0 a])
axis([0 tf -b a])
title('Exact Solution')
hold off

figure(2)
ezsurf(PTb, [0 tf -b 0])
hold on
ezsurf(PTa, [0 tf 0 a])
axis([0 tf -b a])
title('Failed Approximation')
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
zlabel('$T(x,t)$','Interpreter','latex')
colormap gray
hold off
