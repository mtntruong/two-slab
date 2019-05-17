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
clear all; close all; clc;

epsVec = [2 4 6];

Xa = linspace( 0, 1, 20);
Xb = linspace(-1, 0, 20);

Za = zeros(size(Xa));
Zb = zeros(size(Xb));
Ea = zeros(size(Xa));
Eb = zeros(size(Xb));

plot_Za = [];
plot_Zb = [];

for idx = 1 : length(epsVec)
    load(['PT_e' num2str(epsVec(idx))]);

    % Plot at specific 't'

    t0 = 0;
    
	Tb = cos(pi * (x + b)) * cos(pi * y);
    Ta = cos(pi * (x - a)) * cos(pi * y);
    
    for i = 1 : length(Xa)
        Za(i) = subs(PTa,{x, t}, {Xa(i), t0}); % y = 0, t = 0
        Ea(i) = subs( Ta,{x, y}, {Xa(i), 0});  % y = 0
    end
    for i = 1 : length(Xb)
        Zb(i) = subs(PTb,{x, t}, {Xb(i), t0}); % y = 0, t = 0
        Eb(i) = subs( Tb,{x, y}, {Xb(i), 0});  % y = 0
    end
    
    plot_Za(idx,:) = Za;
    plot_Zb(idx,:) = Zb;
end

figure(1)
plot([Xb Xa],[Eb Ea],'k'); hold on;
plot([Xb Xa],[plot_Zb(1,:) plot_Za(1,:)],'ko');
plot([Xb Xa],[plot_Zb(2,:) plot_Za(2,:)],'k+');
plot([Xb Xa],[plot_Zb(3,:) plot_Za(3,:)],'kd');
h = legend('Exact solution',...
       'Material $a$ and $b$ at $10^{-2}$', ...
       'Material $a$ and $b$ at $10^{-4}$', ...
       'Material $a$ and $b$ at $10^{-6}$', ...
       'Location', 'North');
set(h,'Interpreter','latex');
axis([-1 1 -1 1])
xlabel('$x$','Interpreter','latex')
ylabel('$T(x,0)$','Interpreter','latex')
grid on
