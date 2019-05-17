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

function [cA, cB] = generateCMN( indexSet, a, b, c, epsilon, t0, tf, lambda_bar, nu_a, nu_b )

    % Self-explained
    Tb = @(x,y) cos(pi * (x + b)) * cos(pi * y);
    Ta = @(x,y) cos(pi * (x - a)) * cos(pi * y);
    
    % For the ease of typing
    Ym  = @(m,y) sqrt(2) * cos(m * pi * y / c);
    Xan = @(nu,x) cos(nu * (x - a)) / cos(nu * a);
    Xbn = @(nu,x) cos(nu * (x + b)) / cos(nu * b);

    % Grid
	% Partition x-axis by the number of pairs in indexSet to ensure square matrix
    n = size(indexSet,1);
    xla = linspace( 0 + epsilon, a - epsilon, n);
    xlb = linspace(-b + epsilon, 0 - epsilon, n);
    
    % Systems
    LHS_A  = zeros(n, n);
    LHS_B  = zeros(n, n);

    RHS_A  = zeros(n, 1);
    RHS_B  = zeros(n, 1);
    RHS_CA = zeros(n, 1);
    RHS_CB = zeros(n, 1);
    
    % Calculate C_tilde
    % LHS for Ax = B, for each x
    for i = 1 : n
        for j = 1 : n
            % indexSet(j,1), indexSet(j,2) ~ j-th pair m,n
            % nu_a(j), nu_b(j) ~ values of nu corresponding to j-th pair m,n
            LHS_A(i,j) = Xan(nu_a(j), xla(i)) * Ym(indexSet(j,1), 0); % y = 0
            LHS_B(i,j) = Xbn(nu_b(j), xlb(i)) * Ym(indexSet(j,1), 0); % y = 0
        end
    end
    % RHS for calculate C_tilde, at t0
    for i = 1 : n
        RHS_CA(i) = Ta(xla(i), 0); % y = 0
        RHS_CB(i) = Tb(xlb(i), 0); % y = 0
    end
    % Solve for C_tilde at each m
    CmnA = linsolve(LHS_A, RHS_CA);
    CmnB = linsolve(LHS_B, RHS_CB);

    % Calculate T(x,0,tf) and C_mn
    % 
    for i = 1 : n
        % Iterative sum
        RHS_A(i) = 0;
        RHS_B(i) = 0;
        for j = 1 : n
            RHS_A(i) = RHS_A(i) + CmnA(j) * exp(-lambda_bar(j) * tf) * ...
                        Xan(nu_a(j), xla(i)) * Ym(indexSet(j,1), 0); % y = 0
            RHS_B(i) = RHS_B(i) + CmnB(j) * exp(-lambda_bar(j) * tf) * ...
                        Xbn(nu_b(j), xlb(i)) * Ym(indexSet(j,1), 0); % y = 0
        end
        % Add noise
        RHS_A(i) = RHS_A(i) + noise * epsilon;
        RHS_B(i) = RHS_B(i) + noise * epsilon;
    end
    
    % Solve for X = [C_mn], page 7
	cA = linsolve(LHS_A, RHS_A);
	cB = linsolve(LHS_B, RHS_B);
end
