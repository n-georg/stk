% STK_LM_DIFF differentiates a linear model wrt hyper-parameters
%
% CALL: V = stk_lm_diff (LM, X, DIFF)
%
%    computes the partial derivative of the "design matrix" produced
%    the linear model LM at the input points X, with respect to its
%    DIFF-th hyper-parameter.  The result is an N-by-Q, where N is the
%    number of input points, and Q the size of the linear model (number
%    of basis functions).
%
% CALL: [V, LM] = stk_lm_diff (LM, X, DIFF)
%
%    also returns an updated LM object that possibly contains auxiliary
%    results produced during the computation of the gradient.  This
%    provides a mechanism to optimize sequences of computations of
%    partial derivative at a given point.

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

function [V, lm] = stk_lm_diff (lm, x, diff) %#ok<INUSD>

if isa (lm, 'function_handle')
    
    stk_error (['stk_lm_diff cannot be used on linear model objects ' ...
        'represented by function handles, since they have no hyper-' ...
        'parameters'], 'InvalidArgument');
        
else
    
    stk_error (sprintf (['stk_lm_diff is not implemented for models ' ...
        'of class %s.'], class (lm)), 'NotImplemented');
    
end

end % function


%!error stk_lm_diff (@sin, [], 1)
%!error stk_lm_diff (1.2, [], 1)
