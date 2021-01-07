% STK_LM_DIFF [overloas STK function]

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

function [V, lm] = stk_lm_diff (lm, x, diff)  %#ok<INUSD>

if isa (lm, 'stk_lm_')
    
    stk_error (['Classes derived from stk_lm_ must ' ...
        'implement stk_lm_diff.'], 'IncompleteClassImplementation');
    
else
    
    stk_error ('Syntax error', 'SyntaxError');
    
end

end % function
