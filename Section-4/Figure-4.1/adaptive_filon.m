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

function bigA = adaptive_filon(y, a, b, u, g, epsilon)
    
% Nodes and multiplicities
gridSize = 10;
yl = linspace(a, b, gridSize);
ml = 1 * ones(1, gridSize); ml(1) = 2; ml(end) = 2;

% Parameters

omega = g; % g(n)
ell = sum(ml) - 1;

% Symbolic variables

%-- c_q
c = sym(zeros(1, ell));
for i = 0 : ell
    c(i+1) = sym(sprintf('c%d', i));
end

%-- p(y)
p = sym(0);
for i = 0 : ell
    p = p + c(i+1) * y^i;
end

% Solve nonlinear system

%-- Construct matrices
LHS = zeros(ell+1, ell+1);
RHS = zeros(ell+1, 1);

%-- Global counter for linear system indices
sys = 0;

%-- Iteration

% ********** l = 0 **********
%-- main node
%---- p(y) = f(y)
sys = sys + 1;

py = subs(p,y,yl(1));
coeff = coeffExtract(py,c);
LHS(sys,:) = coeff;

u_eps = u + epsilon * noise;
RHS(sys) = eval(subs(u_eps,y,yl(1)));

%---- p'(y) = f'(y) ...
if ml(1) > 1
    temp_P = p;
    temp_U_eps = u + epsilon * noise;
    for j = 1 : (ml(1) - 1)
        temp_P = diff(temp_P);

        sys = sys + 1;

        py = subs(temp_P,y,yl(1));
        coeff = coeffExtract(py,c);
        LHS(sys,:) = coeff;

        RHS(sys) = abs(h_l_omega(epsilon, omega))^-j * delta_j_h0( temp_U_eps, yl(1), j, h_l_omega(epsilon, omega) );
    end
end

% ***** l = 1 ~ nu - 1 ******
if length(yl) > 2
    for i = 2 : (length(yl) - 1)
        %-- main node
        %---- p(y) = f(y)
        sys = sys + 1;
        
        py = subs(p,y,yl(i));
        coeff = coeffExtract(py,c);
        LHS(sys,:) = coeff;

        u_eps = u + epsilon * noise;
        RHS(sys) = eval(subs(u_eps,y,yl(i)));
        
        %---- p'(y) = f'(y) ...
        if ml(i) > 1
            temp_P = p;
            temp_U_eps = u + epsilon * noise;
            for j = 1 : (ml(i) - 1)
                temp_P = diff(temp_P);

                sys = sys + 1;

                py = subs(temp_P,y,yl(i));
                coeff = coeffExtract(py,c);
                LHS(sys,:) = coeff;

                RHS(sys) = abs(h_l_omega(epsilon, omega))^-j * delta_j_hl( temp_U_eps, yl(i), j, h_l_omega(epsilon, omega));
            end
        end
    end
end

% ********* l = nu **********
nu = length(yl);
%-- main node
sys = sys + 1;

py = subs(p,y,yl(nu));
coeff = coeffExtract(py,c);
LHS(sys,:) = coeff;

u_eps = u + epsilon * noise;
RHS(sys) = eval(subs(u_eps,y,yl(nu)));

%---- p'(y) = f'(y) ...
if ml(nu) > 1
    temp_P = p;
    temp_U_eps = u + epsilon * noise;
    for j = 1 : (ml(nu) - 1)
        temp_P = diff(temp_P);

        sys = sys + 1;

        py = subs(temp_P,y,yl(nu));
        coeff = coeffExtract(py,c);
        LHS(sys,:) = coeff;

        RHS(sys) = abs(h_l_omega(epsilon, omega))^-j * nabla_j_hnu( temp_U_eps, yl(nu), j, h_l_omega(epsilon, omega));
    end
end

%-- Display and solve

cq = linsolve(LHS,RHS);

I = zeros(sum(ml),1);
for i = 0 : ell
   I(i+1) = int(y^i * cos(omega * y), a, b);
end

bigA = sum(cq .* I);
