function [yhat predvar] = gpml_gpr_predictor(xstar, covfunc, logtheta, x, y)
% Gaussian process regression, with a named covariance function. It returns
% the prediction and prediction variance of the marginal Gaussian
% predictions. Note that in cases where the covariance function has noise
% contributions, the variance returned in S2 is for noisy test targets; if
% you want the variance of the noise-free latent function, you must
% substract the noise variance.
%
% usage: [yhat predvar] = gpml_gpr_predictor(xstar, covfunc, logtheta, x, y)
%
% where:
%
%   xstar    is a nn by D matrix of test inputs
%   covfunc  is the covariance function
%   logtheta is a (column) vector of log hyperparameters
%   x        is a n by D matrix of training inputs
%   y        is a (column) vector (of size n) of targets
%   yhat     is a (column) vector (of size nn) of prediced means
%   predvar  is a (column) vector (of size nn) of predicted variances
%
% For more help on covariance functions, see "help covFunctions".
%
% (C) copyright 2006 by Carl Edward Rasmussen (2006-03-20).

if ischar(covfunc) % convert to cell if needed
    covfunc = cellstr(covfunc);
end

K = feval(covfunc{:}, logtheta, x); % compute training set covariance matrix
L = chol(K)';                       % cholesky factorization of the covariance
alpha = solve_chol(L',y);

[Kss, Kstar] = feval(covfunc{:}, logtheta, x, xstar); %  test covariances
v = L\Kstar;

yhat    = Kstar' * alpha;
predvar = abs(Kss - sum(v.*v)');

return
