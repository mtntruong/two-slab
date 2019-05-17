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

% PT_e2, PT_e4, and PT_e6 are pre-calculated approximations at different
% epsilons (10^-2, 10^-4, 10^-6, respectively). You can re-produce these
% MAT files by running main.m
clear all; close all; clc

a = 3; b = 5;

xa = linspace(0, a, 20);
xb = linspace(-b, 0, 20);

t0 = 0;

syms x t

for i = 1 : 20
	Ta = eval(getT( 'a', xa(i), a, b ));
	Tb = eval(getT( 'b', xb(i), a, b ));
    line020(i) = eval(subs(subs(Tb, t, t0), x, xb(i)));
    line040(i) = eval(subs(subs(Ta, t, t0), x, xa(i)));
end

figure(1)
plot([xb xa], [line020 line040], 'k');
hold on

load PT_e2
line021 = eval(subs(subs(PTb, t, t0), x, xb));
line041 = eval(subs(subs(PTa, t, t0), x, xa));

plot([xb xa], [line021 line041], 'ko');

load PT_e4
line022 = eval(subs(subs(PTb, t, t0), x, xb));
line042 = eval(subs(subs(PTa, t, t0), x, xa));

plot([xb xa], [line022 line042], 'k+');

load PT_e6
line023 = eval(subs(subs(PTb, t, t0), x, xb));
line043 = eval(subs(subs(PTa, t, t0), x, xa));

plot([xb xa], [line023 line043], 'kd');

L = legend('Exact solution',...
           'Material $a$ and $b$ at $10^{-2}$', ...
           'Material $a$ and $b$ at $10^{-4}$', ...
           'Material $a$ and $b$ at $10^{-6}$', ...
           'Location', 'NorthWest');
set(L,'Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$T(x,0)$','Interpreter','latex')
grid on
