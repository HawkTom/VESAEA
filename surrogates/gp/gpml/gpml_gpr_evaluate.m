function yhat = gpml_gpr_evaluate(xstar, covfunc, logtheta, x, y)
% Gaussian process regression, with a named covariance function. It returns
% the mean of the marginal Gaussian predictions.
%
% usage: yhat = gpml_gpr_evaluate(xstar, covfunc, logtheta, x, y)
%
% where:
%
%   xstar    is the matrix of test inputs
%   covfunc  is the covariance function
%   logtheta is a (column) vector of log hyperparameters
%   x        is the matrix of training inputs
%   y        is a (column) vector (of size n) of targets
%   yhat     is a (column) vector (of size nn) of prediced means
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

yhat = Kstar' * alpha;                                % predicted means

return
