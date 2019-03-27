function [nlml, dnlml] = gpml_gpr_fit(logtheta, covfunc, x, y)
% Gaussian process regression, with a named covariance function. The
% function returns minus the log likelihood and its partial derivatives with
% respect to the hyperparameters. In cases where the covariance function
% has noise contributions, the variance returned in S2 is for noisy test
% targets; if you want the variance of the noise-free latent function, you
% must substract the noise variance.
%
% usage: [nlml dnlml] = gpml_gpr_fit(logtheta, covfunc, x, y)
%
% where:
%
%   logtheta is a (column) vector of log hyperparameters
%   covfunc  is the covariance function
%   x        is a n by D matrix of training inputs
%   y        is a (column) vector (of size n) of targets
%   nlml     is the returned value of the negative log marginal likelihood
%   dnlml    is a (column) vector of partial derivatives of the negative
%            log marginal likelihood wrt each log hyperparameter
%
% For more help on covariance functions, see "help covFunctions".
%
% (C) copyright 2006 by Carl Edward Rasmussen (2006-03-20).

if ischar(covfunc) % convert to cell if needed
    covfunc = cellstr(covfunc);
end

[n, D] = size(x);
if eval(feval(covfunc{:})) ~= size(logtheta, 1)
    error('Error: Number of parameters do not agree with covariance function')
end

K = feval(covfunc{:}, logtheta, x); % compute training set covariance matrix
L = chol(K)';                       % cholesky factorization of the covariance
alpha = solve_chol(L',y);

nlml = 0.5*y'*alpha + sum(log(diag(L))) + 0.5*n*log(2*pi);
dnlml = zeros(size(logtheta));  % set the size of the derivative vector
W = L'\(L\eye(n))-alpha*alpha'; % precompute for convenience
for i = 1:length(dnlml)
    dnlml(i) = sum(sum(W.*feval(covfunc{:}, logtheta, x, i)))/2;
end

return
