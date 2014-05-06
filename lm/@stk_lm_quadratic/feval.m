% FEVAL ... [FIXME: missing documentation]

% Copyright Notice
%
%    Copyright (C) 2012-2014 SUPELEC
%
%    Author:  Julien Bect  <julien.bect@supelec.fr>

% Copying Permission Statement
%
%    This file is part of
%
%            STK: a Small (Matlab/Octave) Toolbox for Kriging
%               (http://sourceforge.net/projects/kriging)
%
%    STK is free software: you can redistribute it and/or modify it under
%    the terms of the GNU General Public License as published by the Free
%    Software Foundation,  either version 3  of the License, or  (at your
%    option) any later version.
%
%    STK is distributed  in the hope that it will  be useful, but WITHOUT
%    ANY WARRANTY;  without even the implied  warranty of MERCHANTABILITY
%    or FITNESS  FOR A  PARTICULAR PURPOSE.  See  the GNU  General Public
%    License for more details.
%
%    You should  have received a copy  of the GNU  General Public License
%    along with STK.  If not, see <http://www.gnu.org/licenses/>.

function z = feval (lm, x) %#ok<INUSL>

x = double (x);

[n, d] = size (x);

z = horzcat (ones (n, 1), x, zeros (n, d * (d + 1) / 2));
k = d + 2;
for i = 1:d
    for j = i:d
        z(:,k) = x(:, i) .* x(:, j);
        k = k + 1;
    end
end

end % function feval


%!test
%! n = 15; d = 4;
%! x = stk_sampling_randunif (n, d);
%! P = feval (stk_lm_quadratic (), x);
%! assert (isequal (size (P), [n, 1 + d * (d + 3) / 2]));
