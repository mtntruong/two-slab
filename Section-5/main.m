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

epsVec = [2 4 6];

for idx = 1 : length(epsVec)

    % Parameters
    
    epsilon = 10^-epsVec(idx);
    
    a = 1; b = 1; c = 1;
    
    kappaA = 0.339;
    kappaB = 0.838;
    
    kA = 1.05;
    kB = 3.42;
    
    t0 = 0; tf = 0.1;
    
    % Initialize
    
    [indexSet, lambda_bar, nu_a, nu_b] = ThetaEps( a, b, c, t0, tf, kA, kB, kappaA, kappaB, epsilon );
    
    [cA, cB] = generateCMN( indexSet, a, b, c, epsilon, t0, tf, lambda_bar, nu_a, nu_b );

    % Symbolic
        
    syms x y t
	PTa = sym(0);
    PTb = sym(0);
    
    % Iteration
    for i = 1 : size(indexSet,1)
        PTa = PTa + cA(i) * exp(lambda_bar(i) * (tf - t)) ...
                * (cos(nu_a(i) * (x - a)) / cos(nu_a(i) * a)) ...
                * (sqrt(2) * cos(indexSet(i,1) * pi * 0)); % y = 0
        PTb = PTb + cB(i) * exp(lambda_bar(i) * (tf - t)) ...
                * (cos(nu_b(i) * (x + b)) / cos(nu_b(i) * b)) ...
                * (sqrt(2) * cos(indexSet(i,1) * pi * 0)); % y = 0
    end
    
    save(['PT_e' num2str(epsVec(idx))]);
    
end
