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

function output = getT( flag, x, a, b )
    
    if (flag == 'a')
        if ( x >= 0 && x <= a/2 )
            output = '-x';
        elseif ( x > a/2 && x <= a )
            output = '-a/2';
        else
            output = 'NULL';
        end
    elseif (flag == 'b')
        if ( x >= -b && x <= -b/2)
            output = '-b/2';
        elseif ( x > -b/2 && x <= 0 )
            output = 'x';
        else
            output = 'NULL';
        end
    end
    
end

