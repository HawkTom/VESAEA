function [srgtSRGT srgtSTT] = srgtsGPFit(srgtOPT)
%Function srgtsGPFit fits the specified Gaussian process model using the
%GPML toolbox of Rasmussen and Williams (2006).
%
%    [srgtSRGT srgtSTT] = srgtsGPFit(srgtOPT)
%
%srgtSRGT is the surrogate structure that contains the following fields:
%* P                    : experimental design matrix.
%* T                    : ouput at P.
%* NbPoints             : number of points in the data set.
%* NbVariables          : number of input variables.
%* GP_CovarianceFunction: covariance function.
%* GP_LogTheta          : covariance function log hyperparameters.
%
%srgtSTT is the state structure that contains the following fields:
%* FIT_Fn    : function handle of the fitting function.
%* FIT_FnVal : value of the loss function (after fitting).
%* NLML      : negative log marginal likelihood.
%* DNLML     : partial derivatives of the negative log marginal likelihood
%            with respect to each log hyperparameter
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
%     [surrogate state] = srgtsGPFit(options)
%
%     surrogate =
%
%                         P: [5x1 double]
%                         T: [5x1 double]
%                  NbPoints: 5
%               NbVariables: 1
%     GP_CovarianceFunction: {'gpml_covSum'  {1x2 cell}}
%               GP_LogTheta: [3x1 double]
%
%     state =
%
%      NLML: 5.9345
%     DNLML: [3x1 double]
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

srgtSRGT.P = srgtOPT.P;
srgtSRGT.T = srgtOPT.T;
[srgtSRGT.NbPoints srgtSRGT.NbVariables] = size(srgtOPT.P);

srgtSRGT.GP_CovarianceFunction = srgtOPT.GP_CovarianceFunction;
srgtSRGT.GP_LogTheta           = srgtOPT.GP_LogTheta0;

switch func2str(srgtOPT.FIT_Fn)
    case 'gpml_minimize'
        % optimization of theta
        srgtSRGT.GP_LogTheta = gpml_minimize(srgtSRGT.GP_LogTheta, 'gpml_gpr_fit', -100, srgtOPT.GP_CovarianceFunction, srgtOPT.P, srgtOPT.T);

        srgtSTT = srgtsFitCreateState(srgtOPT);
        [srgtSTT.NLML srgtSTT.DNLML] = gpml_gpr_fit(srgtSRGT.GP_LogTheta, srgtOPT.GP_CovarianceFunction, srgtOPT.P, srgtOPT.T);
        srgtSTT.FIT_FnVal = srgtSTT.NLML;
        
    case 'gpml_gpr_fit'
        % no optimization for theta
        srgtSTT = srgtsFitCreateState(srgtOPT);
        [srgtSTT.NLML srgtSTT.DNLML] = gpml_gpr_fit(srgtOPT.GP_LogTheta0, srgtOPT.GP_CovarianceFunction, srgtOPT.P, srgtOPT.T);
        srgtSTT.FIT_FnVal = srgtSTT.NLML;
        
    case 'srgtsXVFit'
        [srgtSRGT srgtSTT] = srgtsXVFit(srgtOPT);
        
end

return
