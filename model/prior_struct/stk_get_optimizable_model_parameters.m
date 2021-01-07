% STK_GET_OPTIMIZABLE_PARAMETERS [STK internal]
%
% INTERNAL FUNCTION WARNING:
%
%    This function is currently considered as internal.  STK users that wish to
%    experiment with parameter classes can already overload it, but should be
%    aware that API-breaking changes are likely to happen in future releases.
%
% See also: stk_get_optimizable_parameters

% Copyright Notice
%
%    Copyright (C) 2017, 2018, 2021 CentraleSupelec
%    Copyright (C) 2017 LNE
%
%    Authors:  Remi Stroh   <remi.stroh@lne.fr>
%              Julien Bect  <julien.bect@centralesupelec.fr>

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

function value = stk_get_optimizable_model_parameters (model)

stk_assert_model_struct (model);
model = stk_model_fixlm (model);

% Parameters of the linear model
lmparam = stk_get_optimizable_parameters (model.lm);

% Parameters of the covariance function
covparam = stk_get_optimizable_parameters (model.param);

% Parameters of the noise parameters
noiseparam = stk_get_optimizable_noise_parameters (model);

value = [lmparam; covparam; noiseparam];
   
end % function
