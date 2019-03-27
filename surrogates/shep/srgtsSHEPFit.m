function [srgtSRGT srgtSTT] = srgtsSHEPFit(srgtOPT)
%Function srgtsSHEPFit fits the specified Shepard model.
% 
%    [srgtSRGT srgtSTT] = srgtsSHEPFit(srgtOPT)
% 
%srgtSRGT is the surrogate structure that contains the following fields:
%* P           : experimental design matrix.
%* T           : ouput at P.
%* NbPoints    : number of points in the data set.
%* NbVariables : number of input variables.
%* SHEP_Beta   : coefficients of the local fits.
%* SHEP_Radii  : an array containing the radius of influence for each
%                point.
%
%srgtSTT is the state structure that contains the following fields:
%* SHEP_Error  : 0, if no errors were encountered.
%                1, if N is too small relative to M.
%                2, if any least squares problem is rank deficient.
%                3, if the IRLS subroutine returns an error.
%
%The state structure contains the following field:
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
%     options = srgtsSHEPSetOptions(X, Y);
% 
%     [surrogate state] = srgtsSHEPFit(options)
% 
%     surrogate = 
% 
%               P: [5x1 double]
%               T: [5x1 double]
%        NbPoints: 5
%     NbVariables: 1
%       SHEP_Beta: [-0.6366 -0.6366 -1.1102e-016 0.6366 0.6366]
%      SHEP_Radii: [3.1416 1.5708 1.5708 1.5708 3.1416]
% 
%     state = 
% 
%     SHEP_Error: 0

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

switch srgtOPT.SHEP_Interpolation
    case 'LinearShepard'
        [srgtSRGT.SHEP_Beta, srgtSRGT.SHEP_Radii, srgtSTT.SHEP_Error] = vtechLSHEP(srgtSRGT.NbVariables, srgtSRGT.NbPoints, srgtOPT.P', srgtOPT.T');
        
    case 'Ripple'
        [srgtSRGT.SHEP_Beta, srgtSRGT.SHEP_Radii, srgtSTT.SHEP_Error] = vtechRIPPLE(srgtSRGT.NbVariables, srgtSRGT.NbPoints, srgtOPT.P', srgtOPT.T');
        
end

return
