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

function cf = coeffExtract( symPoly, symArr )
% ** Extract coefficients of symbolic polynomial
% Adding constrains to MATLAB function 'coeffs' to
% make sure it works as intended

% Input:
%   symPoly: symbolic polynomial
%   symArr: list of symbolic variables appear in symPoly

% Output:
%   cf: list of coefficients corresponding to symArr

    cf = zeros(1,length(symArr));
    
    for i = 1 : length(symArr)
        tmp = coeffs(symPoly,symArr(i));
        if length(tmp) > 1
            cf(i) = eval(tmp(2));
        else 
            try
                cf(i) = eval(tmp);
            catch
                cf(i) = 0;
            end
        end
    end
end
