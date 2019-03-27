function predvar = gpml_gpr_predvar(xstar, covfunc, logtheta, x)
% Gaussian process regression, with a named covariance function. It returns
% the variance of the marginal Gaussian predictions. Note that in cases
% where the covariance function has noise contributions, the variance
% returned in S2 is for noisy test targets; if you want the variance of the
% noise-free latent function, you must substract the noise variance.
%
% usage: predvar = gpml_gpr_predvar(xstar, covfunc, logtheta, x)
%
% where:
%
%   logtheta is a (column) vector of log hyperparameters
%   covfunc  is the covariance function
%   x        is a n by D matrix of training inputs
%   y        is a (column) vector (of size n) of targets
%   xstar    is a nn by D matrix of test inputs
%   nlml     is the returned value of the negative log marginal likelihood
%   dnlml    is a (column) vector of partial derivatives of the negative
%            log marginal likelihood wrt each log hyperparameter
%   mu       is a (column) vector (of size nn) of prediced means
%   S2       is a (column) vector (of size nn) of predicted variances
%
% For more help on covariance functions, see "help covFunctions".
%
% (C) copyright 2006 by Carl Edward Rasmussen (2006-03-20).

if ischar(covfunc) % convert to cell if needed
    covfunc = cellstr(covfunc);
end

K = feval(covfunc{:}, logtheta, x); % compute training set covariance matrix
L = chol(K)';                       % cholesky factorization of the covariance

[Kss, Kstar] = feval(covfunc{:}, logtheta, x, xstar);     %  test covariances

v = L\Kstar;
predvar = abs(Kss - sum(v.*v)');

return
