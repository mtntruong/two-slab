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

function generateCN( J, a, b, epsilon, t0, tf, lambda_bar, lambda_a, lambda_b )
    
	syms x y t

    % Grid
    xla = linspace(0 + epsilon, a - epsilon, J);
    xlb = linspace(-b + epsilon, 0 - epsilon, J);
    
    % Systems
    LHS_A  = zeros(J, J);
    LHS_B  = zeros(J, J);
    
    RHS_A  = zeros(J, 1);
    RHS_B  = zeros(J, 1);
    RHS_CA = zeros(J, 1);
    RHS_CB = zeros(J, 1);
    

    % LHS for Ax = B
    for i = 1 : J
        for j = 1 : J
            LHS_A(i,j) = cos(lambda_a(j) * (xla(i) - a)) / cos(lambda_a(j) * a);
            LHS_B(i,j) = cos(lambda_b(j) * (xlb(i) + b)) / cos(lambda_b(j) * b);
        end
    end
    
    % RHS for calculate C_tilde, at t0
    for i = 1 : J
        Ta = eval(getT( 'a', xla(i), a, b ));
        Tb = eval(getT( 'b', xlb(i), a, b ));
        RHS_CA(i) = eval(subs(subs(Ta, t, 0), x, xla(i)));
        RHS_CB(i) = eval(subs(subs(Tb, t, 0), x, xlb(i)));
    end
    
    % Solve for C_tilde
    CnA = linsolve(LHS_A, RHS_CA);
    CnB = linsolve(LHS_B, RHS_CB);
    
    % Calculate T(x,tf)
    for i = 1 : J
        RHS_A(i) = 0;
        RHS_B(i) = 0;
        for n = 1 : J
            RHS_A(i) = RHS_A(i) + CnA(n) * exp(-lambda_bar(n) * tf) * ...
                    cos(lambda_a(n) * (xla(i) - a)) / cos(lambda_a(n) * a);
            RHS_B(i) = RHS_B(i) + CnB(n) * exp(-lambda_bar(n) * tf) * ...
                    cos(lambda_b(n) * (xlb(i) + b)) / cos(lambda_b(n) * b);
        end
        RHS_A(i) = RHS_A(i) + noise * epsilon;
        RHS_B(i) = RHS_B(i) + noise * epsilon;
    end
    
    cA = linsolve(LHS_A, RHS_A);
    cB = linsolve(LHS_B, RHS_B);
    
    save('CA.mat', 'cA');
    save('CB.mat', 'cB');
    
end
