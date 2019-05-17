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

function [indexSet, lbar, nu_a, nu_b] = ThetaEps( a, b, c, t0, tf, kA, kB, kappaA, kappaB, epsilon )

    % Parameters
    beta = 0.05; gamma = 1;

    % Number of the eigenvalues
    Neps = beta * log( epsilon ^ -gamma ) / tf;
    N = floor(Neps);

    lbar = []; nu_a = []; nu_b = [];
    indexSet = []; % indexSet(:,1) ~ m, indexSet(:,2) ~ n
    
    indexSum = -1; % start with m + n = 0, followed by line 16
    while(true)
        flag = false;
        indexSum = indexSum + 1;
        tmpIndex = [];
        for i = 0 : indexSum
            tmpIndex(i+1,1) = i;
            tmpIndex(i+1,2) = indexSum - i;
        end
        
        for i = 1 : size(tmpIndex,1) % iterate through pairs of (m,n)
            % tmpIndex(i,:) ~ (m,n)
            m = tmpIndex(i,1); n = tmpIndex(i,2);
            
            % eigen-function
            F = @(x) kB * x * sin(b*x) * ...
                     cos(a * sqrt(kappaB/kappaA*(x^2) + (kappaB/kappaA - 1) * (m*pi/c)^2)) + ...
                     kA * sqrt(kappaB/kappaA*(x^2) + (kappaB/kappaA - 1) * (m*pi/c)^2) * ...
                     sin(a * sqrt(kappaB/kappaA*(x^2) + (kappaB/kappaA - 1) * (m*pi/c)^2)) * cos(b*x);
                 
            % found_roots = nu_b, find n+1 roots
            tmp_nu_b = find_roots(F, n+1);
            
            % lambda_bar and m,n
            tmp_lbar = kappaB * (tmp_nu_b(n+1)^2 + (m*pi/c)^2); % tmp_nu_b(n+1) ~ nu_bn
            if (tmp_lbar <= N)
                flag = true;
                indexSet = [indexSet; tmpIndex(i,:)];
                tmp_nu_a = sqrt(kappaB/kappaA*(tmp_nu_b(n+1)^2) + (kappaB/kappaA - 1) * (m*pi/c)^2);
                lbar = [lbar tmp_lbar];
                nu_b = [nu_b tmp_nu_b(n+1)];
                nu_a = [nu_a tmp_nu_a];
            end
        end
        
        if (~flag)
            break
        end
    end
end
