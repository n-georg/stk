% STK_PARAM_ESTIM_ [STK internal]
%
% This should become the new stk_param_estim (?).

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

function [model_opt, info] = stk_param_estim_ ...
    (model, xi, zi, param0, lnv0, criterion)

stk_assert_model_struct (model);
% FIXME: Implement stk_param_estim for noise priors too

% Empty is the same as 'not provided'
if nargin < 6
    criterion = [];
    if nargin < 5
        lnv0 = [];
        if nargin < 4
            param0 = [];
        end
    end
end

% Size checking: xi, zi
zi_data = double (zi);
if ~ isequal (size (zi_data), [stk_get_sample_size(xi) 1])
    errmsg = 'zi should be a column, with the same number of rows as xi.';
    stk_error (errmsg, 'IncorrectSize');
end

% Warn about special case: constant response
if (std (zi_data) == 0)
    warning ('STK:stk_param_estim:ConstantResponse', ['Constant-response ' ...
        'data: the output of stk_param_estim is likely to be unreliable.']);
end

% Make sure that lognoisevariance is -inf for noiseless models
if ~ stk_isnoisy (model)
    model.lognoisevariance = -inf;
end

% Default criterion: restricted likelihood (ReML method)
if isempty (criterion)
    criterion = @stk_param_relik;
end

% Make sure that we have a starting point in (param0, lnv0)
[param0, lnv0, do_estim_lnv] = provide_starting_point (model, xi, zi, param0, lnv0);

% Set the starting point
model.param = param0;
model.lognoisevariance = lnv0;

% Get selectors
% (this stuff is ugly, we should simplify by getting rid of do_estim_lnv)
select = stk_param_getblockselectors (model);
lmparam_select = select{1};
covparam_select = select{2};
noiseparam_select = select{3} & do_estim_lnv;
select = [lmparam_select; covparam_select; noiseparam_select];

% NOTE ABOUT SELECTOrS / NaNs
% The situation is currently a little odd, with a focus on estimating
% covariance hyperparameters inherited from STK's history.  More precisely:
%  * the 'linear model' part is not allowed to have (nonlinear) hyperparameters;
%  * the hyperparameters of the covariance function are automatically estimated
%    (all of them), regardless of whether there are NaNs or not in model.param;
%  * the hyperparameters of the noise variance function are estimated (all of
%    them) if there are NaNs in the vector OR if lnv0 is provided.  Note that
%    only the case where lnv is a numerical scalar is currently documented.
%    Also, note that all hyperparameters are estimated, even if only some of
%    them are NaNs.
% FIXME: Clarify this confusing situation...

% Call optimization routine
if nargout >= 2
    [model_opt, info] = stk_param_estim_optim ...
        (model, xi, zi, criterion, select);
else
    model_opt = stk_param_estim_optim ...
        (model, xi, zi, criterion, select);
end

end % function

%#ok<*CTCH,*LERR,*SPWRN,*WNTAG>


function [param0, lnv0, do_estim_lnv] = provide_starting_point ...
    (model, xi, zi, param0, lnv0)

% If param0 is not empty, it means that the user has provided
% a starting point, in which case we must use it.
if ~ isempty (param0)
    
    % Test if param0 contains NaNs of Infs
    assert_is_finite (param0, 'param0');
    % Cast to the appropriate type
    if isnumeric (param0)
        param0 = stk_set_optimizable_parameters (model.param, param0);
    end
    
    % Now take care of lnv0.
    % Same rule: if not empty, we have a user-provided starting point.
    if isempty (lnv0)
        % When lnv0 is not provided, noise variance estimation happens when lnv has NaNs.
        do_estim_lnv = any (isnan (stk_get_optimizable_noise_parameters (model)));
        if do_estim_lnv
            % We have a user-provided starting point for param0 but not for lnv0.
            model.param = param0;
            lnv0 = stk_param_init_lnv (model, xi, zi);
        else
            lnv0 = model.lognoisevariance;
        end
    else
        % When lnv0 is provided, noise variance *must* be estimated
        do_estim_lnv = true;
        % Test if lnv0 contains NaNs of Infs
        assert_is_finite (lnv0, 'lnv0');
        % Cast lnv0 to a variable of the appropriate type (numeric or object)
        lnv0 = stk_set_optimizable_parameters (model.lognoisevariance, lnv0);
    end
    
else  % Otherwise, try stk_param_init to get a starting point
    
    if isempty (lnv0)
        % If needed, stk_param_init will also provide a starting point for lnv0
        % (This is triggered by the presence of NaNs.)
        do_estim_lnv = any (isnan (stk_get_optimizable_noise_parameters (model)));
    else
        % When lnv0 is provided, noise variance *must* be estimated
        do_estim_lnv = true;
        % Test if lnv0 contains NaNs of Infs
        assert_is_finite (lnv0, 'lnv0');        
        % Cast to the appropriate type
        if isnumeric (lnv0)
            lnv0 = stk_set_optimizable_parameters (model.lognoisevariance, param0);
        end
        % Install lnv0 in the model
        model.lognoisevariance = lnv0;
    end
    
    [param0, lnv0] = stk_param_init (model, xi, zi);
    
end % if

end % function


function assert_is_finite (param, param_name)

% Check that there are no NaNs or Infs
param_ = stk_get_optimizable_parameters (param);
if ~ all (isfinite (param_))
    stk_error (sprintf (['Incorrect value for input argument %s.  The components ' ...
        'of the starting point must be neither infinite nor NaN.'], param_name),   ...
        'InvalidArgument');
end

end % function


function select = stk_param_getblockselectors (model)

% FIXME: This function only works for prior model structs.
%        It does not belong here...

stk_assert_model_struct (model);

select = cell (3, 1);

% Linear model parameters
lmparam_size = length (stk_get_optimizable_parameters (model.lm));
select{1} = true (lmparam_size, 1);

% Covariance parameters
covparam_size = length (stk_get_optimizable_parameters (model.param));
select{2} = true (covparam_size, 1);

% Noise parameters
noiseparam_size = length (stk_get_optimizable_noise_parameters (model));
select{3} = true (noiseparam_size, 1);

end % function
