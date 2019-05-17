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

function found_roots = find_roots(F, N)

    interval = [0, 10];
    stop_condition = 0;
    while (~stop_condition)
        M = 1000;
        start_pts = linspace(interval(1), interval(2), M);
        found_roots = [];
        for i = 1 : numel(start_pts)-1
            try
                found_roots(end+1) = fzero(F, [start_pts(i), start_pts(i+1)]);
                if length(found_roots) >= N + 1
                    stop_condition = 1;
                    break
                else
                    interval(2) = interval(2) + 0.3; %d=0.3
                end
            end
        end
    end

end
