% STK_PARAM_PROFLIK computes the profile log-likelihood
%
% CALL: C = stk_param_proflik (MODEL, XI, ZI)
%
%   computes the value C of the opposite of the profile log-likelihood
%   of the MODEL given the data (XI, ZI).
%
% CALL: [C, C_GRAD] = stk_param_proflik (MODEL, XI, ZI)
%
%   also returns the gradient C_GRAD of C with respect to all the
%   optimizable parameters of the model.
%
% REMARK
%
%   This is the "first" profile log-likelihood, obtained by maximizing
%   the log-likelihood with respect to the coefficients of the linear
%   model.  It still depends on the variance parameter, though (if there
%   is one).
%
% EXAMPLE: see paramestim/stk_param_estim.m

% Copyright Notice
%
%    Copyright (C) 2015, 2016, 2018, 2021 CentraleSupelec
%    Copyright (C) 2011-2014 SUPELEC
%
%    Authors:  Julien Bect       <julien.bect@centralesupelec.fr>
%              Emmanuel Vazquez  <emmanuel.vazquez@centralesupelec.fr>

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

function [C, C_grad, C_grad_noiseparam] = stk_param_proflik (model, xi, zi)

zi = double (zi);

% Check the size of zi
n = size (xi, 1);
if ~ isequal (size (zi), [n 1])
    stk_error (['zi must be a column vector, with the same' ...
        'same number of rows as x_obs.'], 'IncorrectSize');
end

LMPRIOR = false;  % Not implemented yet
PARAMPRIOR = isfield (model, 'prior');
NOISEPRIOR = isfield (model, 'noiseprior');

if (nargout >= 2) || LMPRIOR
    % Parameters of the linear model
    lmparam = stk_get_optimizable_parameters (model.lm);
    lmparam_size = length (lmparam);
end

if (nargout >= 2) || PARAMPRIOR
    % Parameters of the covariance function
    covparam = stk_get_optimizable_parameters (model.param);
    covparam_size = length (covparam);
end

if (nargout >= 2) || NOISEPRIOR
    % Parameters of the noise variance function
    [noiseparam, isnoisy] = stk_get_optimizable_noise_parameters (model);
    noiseparam_size = length (noiseparam);
else
    isnoisy = stk_isnoisy (model);
end


%% Compute the (opposite of) the restricted log-likelihood

[K, P] = stk_make_matcov (model, xi);
simple_kriging = (size (P, 2) == 0);

% Choleski factorization: K = U' * U, with upper-triangular U
[U, epsi] = stk_cholcov (K);
UT_TRANSA = struct ('UT', true, 'TRANSA', true);  % Options for linsolve
if (~ isnoisy) && (epsi > 0)
    stk_assert_no_duplicates (xi);
end

% Suffix _tilde for (U')^{-1} * ...
zi_tilde = linsolve (U, zi, UT_TRANSA);

if ~ simple_kriging
    P_tilde = linsolve (U, P, UT_TRANSA);
    beta = (P_tilde' * P_tilde) \ (P_tilde' * zi_tilde);
    % yi = zi - P * beta;
    yi_tilde = zi_tilde - P_tilde * beta;
else
    yi_tilde = zi_tilde;
end

% Compute log (det (K)) using the Cholesky factor
ldetK = 2 * sum (log (diag (U)));

% Compute yi' * K^(-1) * yi
attache = sum (yi_tilde .^ 2);

C = 0.5 * (n * log(2 * pi) + ldetK + attache);


%% Add priors

if PARAMPRIOR
    C = C - stk_distrib_logpdf (model.prior, covparam);
end

if NOISEPRIOR
    C = C - stk_distrib_logpdf (model.noiseprior, noiseparam);
end


%% Compute gradient

if nargout >= 2
    
    % Compute K^(-1) from the Cholesky factor
    if exist ('OCTAVE_VERSION', 'builtin') == 5
        K_inv = chol2inv (U);  % = K^(-1)
    else
        % Matlab does not have chol2inv
        % TODO: Write a mex to call LAPACK's dpotri
        F = linsolve (U, eye (n), UT_TRANSA);
        K_inv = F' * F;  % = K^(-1)
    end
    
    w = linsolve (U, yi_tilde, struct ('UT', true));
    
    % Gradient wrt parameters of the linear model
    C_grad_lmparam = zeros (lmparam_size, 1);
    if lmparam_size > 0
        lm_ = model.lm;
        for diff = 1:lmparam_size
            [V, lm_] = stk_lm_diff (lm_, xi, diff);
            C_grad_lmparam(diff) = - w' * V * beta;
        end
        if LMPRIOR
            % FIXME: Implement!
            stk_error ('LMPRIOR is not implemented yet.', ...
                'NotImplemented');  %#ok<UNRCH>
        end
    end
    
    % Gradient wrt parameters of the covariance function
    C_grad_covparam = zeros (covparam_size, 1);
    if covparam_size > 0
        for diff = 1:covparam_size
            V = feval (model.covariance_type, model.param, xi, xi, diff);
            C_grad_covparam(diff) = 1/2 * (sum (sum (K_inv .* V)) - w' * V * w);
        end
        if PARAMPRIOR
            C_grad_covparam = C_grad_covparam ...
                - stk_distrib_logpdf_grad (model.prior, covparam);
        end
    end
    
    % Gradient wrt parameters of the noise model
    C_grad_noiseparam = zeros (noiseparam_size, 1);
    if noiseparam_size > 0
        for diff = 1:noiseparam_size
            V = stk_covmat_noise (model, xi, [], diff);
            C_grad_noiseparam(diff) = 1/2 * (sum (sum (K_inv .* V)) - w' * V * w);
        end
        if NOISEPRIOR
            C_grad_noiseparam = C_grad_noiseparam ...
                - stk_distrib_logpdf_grad (model.noiseprior, noiseparam);
        end
    end
    
end


%% Construct full gradient (if needed)

switch nargout
    
    case {0, 1}
        % Nothing to do, the gradient was not requested
        
    case 2
        % Recommended syntax, with the full gradient as the second output
        C_grad = [C_grad_lmparam; C_grad_covparam; C_grad_noiseparam];
        
    case 3
        % Old (deprecated) syntax, from a time where linear models were
        % not allowed to have additional parameters
        assert (isempty (C_grad_lmparam));
        C_grad = C_grad_covparam;
        
    otherwise
        stk_error ('Too many input arguments.', 'TooManyInputArgs');
        
end % switch

end % function


%!shared f, xi, zi, NI, model, C, dC1, dC2
%!
%! f = @(x)(- (0.8 * x(:, 1) + sin (5 * x(:, 2) + 1) ...
%!          + 0.1 * sin (10 * x(:, 3))));
%! DIM = 3;  NI = 20;  box = repmat ([-1.0; 1.0], 1, DIM);
%! xi = stk_sampling_halton_rr2 (NI, DIM, box);
%! zi = stk_feval (f, xi);
%!
%! SIGMA2 = 1.0;  % variance parameter
%! NU     = 4.0;  % regularity parameter
%! RHO1   = 0.4;  % scale (range) parameter
%!
%! model = stk_model('stk_materncov_aniso');
%! model.param = log([SIGMA2; NU; 1/RHO1 * ones(DIM, 1)]);

%!error [C, dC1, dC2] = stk_param_proflik ();
%!error [C, dC1, dC2] = stk_param_proflik (model);
%!error [C, dC1, dC2] = stk_param_proflik (model, xi);
%!test  [C, dC1, dC2] = stk_param_proflik (model, xi, zi);

% %!test
% %! TOL_REL = 0.01;
% %! assert (stk_isequal_tolrel (C, 21.6, TOL_REL));
% %! assert (stk_isequal_tolrel (dC1, [4.387 -0.1803 0.7917 0.1392 2.580]', TOL_REL));
% %! assert (isequal (dC2, zeros (0, 1)));

% %!shared xi, zi, model, TOL_REL
% %! xi = [-1 -.6 -.2 .2 .6 1]';
% %! zi = [-0.11 1.30 0.23 -1.14 0.36 -0.37]';
% %! model = stk_model ('stk_materncov_iso');
% %! model.param = log ([1.0 4.0 2.5]);
% %! model.lognoisevariance = log (0.01);
% %! TOL_REL = 0.01;

% %!test  % Another simple 1D check
% %! [C, dC1, dC2] = stk_param_relik (model, xi, zi);
% %! assert (stk_isequal_tolrel (C, 6.327, TOL_REL));
% %! assert (stk_isequal_tolrel (dC1, [0.268 0.0149 -0.636]', TOL_REL));
% %! assert (stk_isequal_tolrel (dC2, -1.581e-04, TOL_REL));

% %!test  % Same 1D test with simple kriging
% %! model.lm = stk_lm_null;
% %! [C, dC1, dC2] = stk_param_relik (model, xi, zi);
% %! assert (stk_isequal_tolrel (C, 7.475, TOL_REL));
% %! assert (stk_isequal_tolrel (dC1, [0.765 0.0238 -1.019]', TOL_REL));
% %! assert (stk_isequal_tolrel (dC2, 3.0517e-03, TOL_REL));

%!test  % Check the gradient on a 2D test case
%!
%! f = @stk_testfun_braninhoo;
%! DIM = 2;
%! BOX = [[-5; 10], [0; 15]];
%! NI = 20;
%! TOL_REL = 1e-2;
%! DELTA = 1e-6;
%!
%! model = stk_model ('stk_materncov52_iso', DIM);
%! model.param = [1 1];
%!
%! xi = stk_sampling_halton_rr2 (NI, DIM, BOX);
%! zi = stk_feval (f, xi);
%!
%! for range = [0.3 2 10]
%!     model.param(2) = - log (range);
%!     for diff = 1:2
%!         assert (stk_test_critgrad ...
%!             (@stk_param_proflik, model, xi, zi, diff, 1e-6));
%!     end
%! end
