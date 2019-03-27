function srgtOPT = srgtsGPSetOptions(P, T, FIT_Fn, ...
    GP_CovarianceFunction, GP_LogTheta0)
%Function srgtsGPSetOptions creates the SURROGATES Toolbox option
%structure for Gaussian process models. This structure contains the
%following fiels:
%
%* GENERAL PARAMETERS
%
%   SRGT   - Identifier of the surrogate technique: 'GP'.
%   P      - NPOINTS-by-NDV matrix, where NPOINTS is the number of points
%            of the sample and NDV is the number of design variables.
%            Default: Empty matrix.
%   T      - NPOINTS-by-1 vector of responses on the P matrix points.
%            Default: Empty vector.
%   FIT_Fn - Function handle of the fitting function (which is used to
%            optimize GP_LogTheta). [@gpml_minimize | @gpml_gpr_fit |
%            @srgtsXVFit]. Default: @gpml_gpr_fit.
%
%* GAUSSIAN PROCESS PARAMETERS
%
%   GP_CovarianceFunction - covariance function. [ string | 'gpml_covConst'
%                           | 'gpml_covLINard' | 'gpml_covLINone' |
%                           'gpml_covMatern3iso' | 'gpml_covMatern5iso' |
%                           'gpml_covNNone' | 'gpml_covNoise' |
%                           'gpml_covPeriodic' | 'gpml_covProd' |
%                           'gpml_covRQard' | 'gpml_covRQiso' |
%                           'gpml_covSEard' | 'gpml_covSEiso' |
%                           'gpml_covSum']. Default: 'gpml_covSEard'.
%   GP_LogTheta0          - Initial guess for LogTheta (covariance function
%                           log hyperparameters). Default:
%                           [log(NPOINTS^(-1/NDV))*ones(1, NDV)/10; 0].
%   GP_LowerBound         - Lower bound for LogTheta. Default: Empty vector.
%   GP_UpperBound         - Upper bound for LogTheta. Default: Empty vector.
%
%The SURROGATES Toolbox uses the GPML toolbox of Rasmussen and Williams
%(2006) to execute the Gaussian process algorithm.
% 
%This is how you can use srgtsGPSetOptions:
% 
%     OPTIONS = srgtsGPSetOptions: creates a structure with the empty
%     parameters.
%
%     OPTIONS = srgtsGPSetOptions(P, T): Given the sampled data P (input
%     variables) and T (output variables), it creates a structure with
%     default parameters used for all not specified fields.
%
%     OPTIONS = srgtsGPSetOptions(P, T, FIT_Fn, GP_CovarianceFunction, ...
%     GP_LogTheta0, GP_LowerBound, GP_UpperBound):  creates a structure
%     with each of the specified fields.
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
%     options = srgtsGPSetOptions(X, Y)
%
%     options =
%
%                      SRGT: 'GP'
%                         P: [5x1 double]
%                         T: [5x1 double]
%     GP_CovarianceFunction: 'gpml_covSEiso'
%              GP_LogTheta0: [3x1 double]
%          GP_OptimizeTheta: 0
%
%REFERENCES:
%
%Rasmussen CE and Williams CKI, Gaussian Processes for Machine Learning,
%The MIT Press, 2006.
%Available at: http://www.gaussianprocess.org/gpml/.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
srgtOPT.SRGT = 'GP';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options
switch nargin
    case 0
        srgtOPT.P = [];
        srgtOPT.T = [];
        
        srgtOPT.FIT_Fn = [];
        
        srgtOPT.GP_CovarianceFunction = [];
        srgtOPT.GP_LogTheta0          = [];
        srgtOPT.GP_LowerBound       = [];
        srgtOPT.GP_UpperBound       = [];
    case 2
        [npoints nvariables] = size(P);

        srgtOPT.P = P;
        srgtOPT.T = T;
        
        srgtOPT.FIT_Fn = @gpml_gpr_fit;
        
        srgtOPT.GP_CovarianceFunction = 'gpml_covSEard';
        srgtOPT.GP_LogTheta0          = [log(npoints^(-1/nvariables))*ones(nvariables, 1)/10; 0];
        srgtOPT.GP_LowerBound       = [];
        srgtOPT.GP_UpperBound       = [];
    otherwise
        srgtOPT.P = P;
        srgtOPT.T = T;

        srgtOPT.FIT_Fn = FIT_Fn;
        
        srgtOPT.GP_CovarianceFunction = GP_CovarianceFunction;
        srgtOPT.GP_LogTheta0          = GP_LogTheta0;
        srgtOPT.GP_LowerBound       = [];
        srgtOPT.GP_UpperBound       = [];
end

return
