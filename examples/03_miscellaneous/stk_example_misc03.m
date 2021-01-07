% STK_EXAMPLE_MISC03  How to deal with a simple seasonality effect

% Copyright Notice
%
%    Copyright (C) 2016, 2021 CentraleSupelec
%    Copyright (C) 2014 SUPELEC
%
%    Author:  Julien Bect  <julien.bect@centralesupelec.fr>

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

stk_disp_examplewelcome


%% Sinusoid + noise

t_obs = (0:0.05:12)';
S_obs = sin (t_obs + 0.3) + 0.1 * randn (size (t_obs));

stk_figure ('stk_example_misc03');
plot (t_obs, S_obs, 'k.', 'DisplayName', 'data');
stk_labels ('month number t', 'sunspots S');
stk_legend ();

t = (0:0.01:30)';


%% Gaussian process model with constant prior mean

model = stk_model (@stk_materncov52_iso);
model.lognoisevariance = nan;

% Initial guess for the parameters of the Matern covariance
[param0, lnv0] = stk_param_init (model, t_obs, S_obs);

% Estimate the parameters
[model.param, model.lognoisevariance] = stk_param_estim ...
    (model, t_obs, S_obs, param0, lnv0);

% Carry out the kriging prediction
S_post = stk_predict (model, t_obs, S_obs, t);

% Display the result
hold on;  plot (t, S_post.mean, 'r-', 'DisplayName', 'constant mean');


%% Gaussian process model with (known) seasonality

% Periodicity assumed to be known
T0 = 2 * pi;

% Construct a prior model with sinusoidal trend
model2 = stk_model (@stk_materncov52_iso);
model2.lm = @(t)([ones(length(t),1) sin(2*pi*t/T0) cos(2*pi*t/T0)]);
model2.lognoisevariance = nan;

% Initial guess for the parameters of the Matern covariance
[param0, lnv0] = stk_param_init (model2, t_obs, S_obs);

% Estimate the parameters
[model2.param, model2.lognoisevariance] = ...
    stk_param_estim (model2, t_obs, S_obs, param0, lnv0);

% Carry out the kriging prediction
S_post = stk_predict (model2, t_obs, S_obs, t);

% Display the result
hold on;  plot (t, S_post.mean, 'g-', ...
    'LineWidth', 3, 'DisplayName', 'known \omega');


%% Gaussian process model with (partially unknown) seasonality

% We assume that the period is unknown, but the number of harmonics
% (equal to one) is known.

% Initial guess and search interval for omega
omega0 = 1.2;  % The true value is omega = 1.0
omega_lims = [0.1; 5.0];

% Construct a prior model with sinusoidal trend
model3 = stk_model ('stk_materncov52_iso');
model3.lm = stk_lm_sincos (omega0, omega_lims);
model3.lognoisevariance = nan;

% FIXME: Make it possible to include a constant term as well
%        (sum of linear model objects !)

% Initial guess for the parameters of the Matern covariance
[param0, lnv0] = stk_param_init (model3, t_obs, S_obs);

% NOTE: stk_param_init completely ignore the existing lm field
%       and uses an stk_lm_constant term instead (FIXME ?)

% FIXME: "old" [param, lnv] syntax, cf. stk_param_estim

% Estimate the parameters  (NEW EXPERIMENTAL SYNTAX !!!)
model3 = stk_param_estim_ (model3, ...
    t_obs, S_obs, param0, lnv0, @stk_param_proflik);

% Carry out the kriging prediction
S_post = stk_predict (model3, t_obs, S_obs, t);

% Display the result
hold on;  plot (t, S_post.mean, 'b--', ...
    'LineWidth', 2, 'DisplayName', 'unknown \omega');

stk_legend ();


%% Display models

model
model2
model3

%#ok<*NOPTS>

%!test stk_example_misc03;  close all;
