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

function generateCN( J, a, b, C, epsilon, t0, tf, lambda_a, lambda_b )
    
    % Grid
    
    xla = linspace(0 + epsilon, a - epsilon, J);
    xlb = linspace(-b + epsilon, 0 - epsilon, J);
    
    % Symbolic
    
    syms x y t
    
    Ta = C * cos(lambda_a(2) * (x - a)) * exp(-t) / ...
        ((a + b) * cos (lambda_a(2) * a));
    Tb = C * cos(lambda_b(2) * (x + b)) * exp(-t) / ...
        ((a + b) * cos (lambda_b(2) * b));
    
    LHS_A = zeros(J, J);
    LHS_B = zeros(J, J);
    
    RHS_A = eval(subs(subs(Ta + epsilon * noise, t, tf), x, xla))';
    RHS_B = eval(subs(subs(Tb + epsilon * noise, t, tf), x, xlb))';
    
    for i = 1 : J
        for j = 1 : J
            LHS_A(i,j) = cos(lambda_a(j) * (xla(i) - a)) / cos(lambda_a(j) * a);
            LHS_B(i,j) = cos(lambda_b(j) * (xlb(i) + b)) / cos(lambda_b(j) * b);
        end
    end
    
    cA = linsolve(LHS_A, RHS_A);
    cB = linsolve(LHS_B, RHS_B);
    
    save('CA.mat', 'cA');
    save('CB.mat', 'cB');
    
end
