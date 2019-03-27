function srgtOPT = srgtsSHEPSetOptions(P, T, SHEP_Interpolation)
%Function srgtsSHEPSetOptions creates the SURROGATES Toolbox option
%structure for Shepard models. This structure contains
%the following fields:
%
%* GENERAL PARAMETERS
%
%   SRGT - Identifier of the surrogate technique: 'SHEP'.
%   P    - NPOINTS-by-NDV matrix, where NPOINTS is the number of points of
%          the sample and NDV is the number of design variables.
%          Default: Empty matrix.
%   T    - NPOINTS-by-1 vector of responses on the P matrix points.
%          Default: Empty vector.
%
%* SHEPARD PARAMETERS
%   SHEP_Interpolation - Interpolation scheme [ string | 'LinearShepard' |
%                        'Ripple' ]. Default: 'LinearShepard'
%
%This is how you can use srgtsSHEPSetOptions:
% 
%     OPTIONS = srgtsSHEPSetOptions: creates a structure with the empty
%     parameters.
%
%     OPTIONS = srgtsSHEPSetOptions(P, T): Given the sampled data P (input
%     variables) and T (output variables), it creates a structure with
%     default parameters used for all not specified fields.
%
%     OPTIONS = srgtsSHEPSetOptions(P, T, SHEP_Interpolation): it creates a
%     structure with each of the specified fields.
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
%     options = srgtsSHEPSetOptions(X, Y)
%
%     options =
%
%                   SRGT: 'SHEP'
%                      P: [5x1 double]
%                      T: [5x1 double]
%     SHEP_Interpolation: 'LinearShepard'

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
srgtOPT.SRGT = 'SHEP';

srgtOPT.P = P;
srgtOPT.T = T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options
if nargin == 2
    srgtOPT.SHEP_Interpolation = 'LinearShepard';
else
    srgtOPT.SHEP_Interpolation = SHEP_Interpolation;
end

return
