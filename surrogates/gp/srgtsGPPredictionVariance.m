function predvar = srgtsGPPredictionVariance(x, srgtSRGT)
%Function srgtsGPPredictionVariance returns the estimated prediction
%variance of a Gaussian process model. Thus, for example:
%
%     PREDVAR = srgtsGPPredictionVariance(X, SURROGATE): returns the
%     prediction variance PREDVAR estimated by the Gaussian process model SURROGATE
%     at all X sites. X can be either a single row vector (single point,
%     with each column representing a variable) or a matrix (each row
%     represents a point). PREDVAR is an NPOINTS-by-1 vector.
%
%Example:
%     % basic information about the problem
%     myFN = @cos;  % this could be any user-defined function
%     designspace = [0;     % lower bound
%                    2*pi]; % upper bound
%
%     % create DOE
%     npoints = 5;
%     X = linspace(designspace(1), designspace(2), npoints)';
%
%     % evaluate analysis function at X points
%     Y = feval(myFN, X);
%
%     % fit surrogate models
%     options = srgtsGPSetOptions(X, Y);
% 
%     [surrogate state] = srgtsGPFit(options);
%
%     % create test points
%     Xtest = linspace(designspace(1), designspace(2), 100)';
%
%     % evaluate surrogate at Xtest
%     PredVar = srgtsGPPredictionVariance(Xtest, surrogate);
%
%     plot(X, zeros(npoints, 1), 'o', ...
%          Xtest, PredVar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

predvar = gpml_gpr_predvar(x, srgtSRGT.GP_CovarianceFunction, srgtSRGT.GP_LogTheta, srgtSRGT.P);

return
