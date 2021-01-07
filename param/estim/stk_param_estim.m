% STK_PARAM_ESTIM estimates the parameters of a covariance function
%
% CALL: PARAM = stk_param_estim (MODEL, XI, YI, PARAM0)
% CALL: [PARAM, LNV] = stk_param_estim (MODEL, XI, YI, PARAM0)
%
%   estimates the parameters PARAM of the covariance function in MODEL
%   from the data (XI, YI) using the restricted maximum likelihood (ReML)
%   method.  The value PARAM0 is used as a starting point for local
%   optimization.
%
%   The observations are assumed to be noisy if MODEL.lognoisevariance is
%   not -inf. In this case, the variance of the noise is estimated if
%   MODEL.lognoisevariance is nan, and assumed known otherwise. The
%   estimated log-variance is returned as the second output argument LNV
%   (equal to MODEL.lognoisevariance when it is assumed to be known).
%
% CALL: PARAM = stk_param_estim (MODEL, XI, YI)
% CALL: [PARAM, LNV] = stk_param_estim (MODEL, XI, YI)
%
%   does the same thing but uses stk_param_init to provide a starting value
%   automatically.
%
% CALL: [PARAM, LNV] = stk_param_estim (MODEL, XI, YI, PARAM0, LNV0)
%
%   additionally provides an initial guess LNV0 for the logarithm of the
%   noise variance. In this case the observations are automatically assumed
%   to be noisy, and the value of MODEL.lognoisevariance is ignored.
%
% CALL: PARAM = stk_param_estim (MODEL, XI, YI, PARAM0, [], CRIT)
% CALL: [PARAM, LNV] = stk_param_estim (MODEL, XI, YI, PARAM0, LNV0, CRIT)
%
%   uses the estimation criterion CRIT instead of the default ReML criterion.
%
% EXAMPLES: see, e.g., stk_example_kb02, stk_example_kb03, stk_example_kb04,
%           stk_example_kb06, stk_example_misc02

% Copyright Notice
%
%    Copyright (C) 2015-2018, 2020, 2021 CentraleSupelec
%    Copyright (C) 2017 LNE
%    Copyright (C) 2014 A. Ravisankar
%    Copyright (C) 2011-2014 SUPELEC
%
%    Authors:  Julien Bect        <julien.bect@centralesupelec.fr>
%              Emmanuel Vazquez   <emmanuel.vazquez@centralesupelec.fr>
%              Ashwin Ravisankar  <ashwinr1993@gmail.com>
%              Remi Stroh         <remi.stroh@lne.fr>

% Copying Permission Statement
%
%    This file is part of
%
%            STK: a Small (Matlab/Octave) Toolbox for Kriging
%               (https://github.com/stk-kriging/stk/)
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

function [param_opt, lnv_opt, info] = stk_param_estim (varargin)

% TODO (?): Remove this function entirely, use stk_param_estim_ instead,
%           and jump to major release 3.x

switch nargout
    
    case {0, 1}
        model_opt = stk_param_estim_ (varargin{:});
        param_opt = model_opt.param;
        
    case 2
        model_opt = stk_param_estim_ (varargin{:});
        param_opt = model_opt.param;
        lnv_opt = model_opt.lognoisevariance;
        
    case 3
        [model_opt, info] = stk_param_estim_ (varargin{:});
        param_opt = model_opt.param;
        lnv_opt = model_opt.lognoisevariance;
        
    otherwise
        stk_error ('Too many input arguments.', 'TooManyInputArgs');
        
end

end % function

%#ok<*CTCH,*LERR,*SPWRN,*WNTAG>


%!shared f, xi, zi, NI, param0, param1, model
%!
%! f = @(x)(- (0.8 * x + sin (5 * x + 1) + 0.1 * sin (10 * x)) );
%! DIM = 1;  NI = 20;  box = [-1.0; 1.0];
%! xi = stk_sampling_regulargrid (NI, DIM, box);
%!
%! SIGMA2 = 1.0;  % variance parameter
%! NU     = 4.0;  % regularity parameter
%! RHO1   = 0.4;  % scale (range) parameter
%! param0 = log ([SIGMA2; NU; 1/RHO1]);
%!
%! model = stk_model (@stk_materncov_iso);

%!test  % noiseless
%! zi = stk_feval (f, xi);
%! param1 = stk_param_estim (model, xi, zi, param0);
%! assert (isequal (size (param1), size (param0)))

% We cannot assume a DETERMINISTIC optimization algorithm
% (for some reason, Octave's sqp is not exactly deterministic)

%!test  % same thing, with empty lnv0 (ok)
%! param2 = stk_param_estim (model, xi, zi, param0, []);
%! assert (stk_isequal_tolrel (param2, param1, 1e-2))

%!error  % same thing, with lnv0 == NaN (not ok as a starting point)
%! param2 = stk_param_estim (model, xi, zi, param0, nan);

%!error  % same thing, with lnv0 == -inf (not ok as a starting point)
%! param2 = stk_param_estim (model, xi, zi, param0, -inf);

%!test  % same thing, with explicit value for 'criterion'
%! param2 = stk_param_estim (model, xi, zi, param0, [], @stk_param_relik);
%! assert (stk_isequal_tolrel (param1, param2, 1e-2))

%!test  % noiseless
%! zi = stk_feval (f, xi);
%! param1 = stk_param_estim (model, xi, zi, param0);
%! param2 = stk_param_estim (model, xi, zi, param0, [], @stk_param_relik);
%! % We cannot assume a DETERMINISTIC optimization algorithm
%! % (for some reason, Octave's sqp is not exactly deterministic)
%! assert (stk_isequal_tolrel (param1, param2, 1e-2))

%!test  % noisy
%! NOISE_STD_TRUE = 0.1;
%! NOISE_STD_INIT = 1e-5;
%! zi = zi + NOISE_STD_TRUE * randn(NI, 1);
%! model.lognoisevariance = 2 * log(NOISE_STD_INIT);
%! [param, lnv] = stk_param_estim ...
%!    (model, xi, zi, param0, model.lognoisevariance);

% Incorrect number of input arguments
%!error param = stk_param_estim ()
%!error param = stk_param_estim (model);
%!error param = stk_param_estim (model, xi);

%!test % Constant response
%! model = stk_model (@stk_materncov52_iso);
%! n = 10;  x = stk_sampling_regulargrid (n, 1, [0; 1]);  z = ones (size (x));
%! param = stk_param_estim (model, x, z);
%! assert ((all (isfinite (param))) && (length (param) == 2));
