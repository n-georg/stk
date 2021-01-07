% STK_LM_SINCOS creates a 'sin+cos' linear model object
%
% TODO: document me!

% Copyright Notice
%
%    Copyright (C) 2021 CentraleSupelec
%
%    Author:  Julien Bect  <julien.bect@centralesupelec.fr>

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

function lm = stk_lm_sincos (omega, omega_lims)

if nargin == 0
    
    lm = struct (                    ...
        'angular_frequency',  1.0,   ...
        'param',              0.0,   ...
        'param_min',          -inf,  ...
        'param_max',          +inf   );
    
else
    
    % FIXME: Check input size and type
    
    assert (omega_lims(1) <= omega_lims(2));
    
    lm = struct (                                   ...
        'angular_frequency',  omega,                ...
        'param',              log (omega),          ...
        'param_min',          log (omega_lims(1)),  ...
        'param_max',          log (omega_lims(2))   );
    
end

lm = class (lm, 'stk_lm_sincos', stk_lm_ ());

end  % function


%!test stk_test_class ('stk_lm_sincos')
