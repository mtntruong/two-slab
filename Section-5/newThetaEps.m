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

function [lbar, nu_a, nu_b] = newThetaEps( M, N, a, b, c, t0, tf, kA, kB, kappaA, kappaB, epsilon )

    % Parameters
    beta = 0.01; gamma = 1;

    lbar = []; nu_a = []; nu_b = [];
        
    for i = 0 : M % iterate through m
        m = i;
        
        % eigen-function
        F = @(x) kB * x * sin(b*x) * ...
                 cos(a * sqrt(kappaB/kappaA*(x^2) + (kappaB/kappaA - 1) * (m*pi/c)^2)) + ...
                 kA * sqrt(kappaB/kappaA*(x^2) + (kappaB/kappaA - 1) * (m*pi/c)^2) * ...
                 sin(a * sqrt(kappaB/kappaA*(x^2) + (kappaB/kappaA - 1) * (m*pi/c)^2)) * cos(b*x);
        tmp_nu_b = find_roots(F, N+1);
        for j = 0 : N
            n = j;
            tmp_lbar = kappaB * (tmp_nu_b(n+1)^2 + (m*pi/c)^2); % tmp_nu_b(n+1) ~ nu_bn
            tmp_nu_a = sqrt(kappaB/kappaA*(tmp_nu_b(n+1)^2) + (kappaB/kappaA - 1) * (m*pi/c)^2);
            lbar(m+1,n+1) = tmp_lbar;
            nu_b(m+1,n+1) = tmp_nu_b(n+1);
            nu_a(m+1,n+1) = tmp_nu_a;
        end
    end

end
