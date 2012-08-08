% STK_SF_MATERN32 computes the Matern correlation function of order 3/2.
%
% CALL: K = stk_sf_matern32(H)
%
%    computes the value of the Matern correlation function of order 3/2 at
%    distance H. Note that the Matern correlation function is a valid
%    correlation function for all dimensions.
%
% CALL: K = stk_sf_matern32(H, DIFF)
%    
%    computes the derivative of the Matern correlation function of order 3/2, at
%    distance H, with respect the distance H if DIFF is equal to 1. (If DIFF is 
%    equal to -1, this is the same as K = stk_sf_matern32(H).)
%
% See also: stk_sf_matern, stk_sf_matern52

%                  Small (Matlab/Octave) Toolbox for Kriging
%
% Copyright Notice
%
%    Copyright (C) 2011, 2012 SUPELEC
%    Version:   1.1
%    Authors:   Julien Bect       <julien.bect@supelec.fr>
%               Emmanuel Vazquez  <emmanuel.vazquez@supelec.fr>
%    URL:       http://sourceforge.net/projects/kriging/
%
% Copying Permission Statement
%
%    This  file is  part  of  STK: a  Small  (Matlab/Octave) Toolbox  for
%    Kriging.
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
%
function k = stk_sf_matern32(h, diff)

stk_narginchk(1, 2);

% default: compute the value (not a derivative)
if (nargin<2), diff = -1; end

Nu = 3/2;
C  = 2 * sqrt(Nu);   % dt/dh
t  = C * abs(h);

switch diff,
    
    case -1, % value of the covariance function
        
        k = (1 + t) .* exp(-t);
        
    case 1, % derivative wrt h
            
        k = - C * t .* exp(-t);

    otherwise
        
        error('incorrect value for diff.');
        
end

end


%%%%%%%%%%%%%
%%% tests %%%
%%%%%%%%%%%%%

%!shared h, diff
%! h = 1.0; diff = -1;

%!error stk_sf_matern32();
%!test  stk_sf_matern32(h);
%!test  stk_sf_matern32(h, diff);
%!error stk_sf_matern32(h, diff, pi);

%!test %% h = 0.0 => correlation = 1.0
%! x = stk_sf_matern32(0.0);
%! assert(stk_isequal_tolrel(x, 1.0, 1e-8));

%!test %% consistency with stk_sf_matern: function values
%! for h = 0.1:0.1:2.0,
%!   x = stk_sf_matern(3/2, h);
%!   y = stk_sf_matern32(h);
%!   assert(stk_isequal_tolrel(x, y, 1e-8));
%! end

%!test %% consistency with stk_sf_matern: derivatives
%! for h = 0.1:0.1:2.0,
%!   x = stk_sf_matern(3/2, h, 2);
%!   y = stk_sf_matern32(h, 1);
%!   assert(stk_isequal_tolrel(x, y, 1e-8));
%! end