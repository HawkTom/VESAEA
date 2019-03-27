function CMatrix = srgtsWASComputeCMatrix(P, ErrorMatrix, flag)
%Function srgtsWASComputeCMatrix computes the C matrix used in some
%weighted average surrogate schemes. Thus, for example:
%
%     CMATRIX = srgtsWASComputeCMatrix(P, ERRORMATRIX, FLAG):
%     computes the C matrix CMATRIX based on the given P data points and
%     the respective ERRORMATRIX error matrix. It is important to say that
%     all P(i) must range from 0 to 1. FLAG controls if it is necessary to
%     correct (FLAG = 1) or not (FLAG = 0) ERRORMATRIX to take into account
%     the effect of trapezoidal integration.
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
%     % kriging
%     srgtOPTKRG  = srgtsKRGSetOptions(X, Y);
%     srgtSRGTKRG = srgtsKRGFit(srgtOPTKRG);
%     [PRESSRMS_KRG, eXV_KRG] = srgtsCrossValidation(srgtOPTKRG);
% 
%     % support vector regression
%     srgtOPTSVR  = srgtsSVRSetOptions(X, Y);
%     srgtSRGTSVR = srgtsSVRFit(srgtOPTSVR);
%     [PRESSRMS_SVR, eXV_SVR] = srgtsCrossValidation(srgtOPTSVR);
% 
%     % shepard
%     srgtOPTSHEP  = srgtsSHEPSetOptions(X, Y);
%     srgtSRGTSHEP = srgtsSHEPFit(srgtOPTSHEP);
%     [PRESSRMS_SHEP, eXV_SHEP] = srgtsCrossValidation(srgtOPTSHEP);
%
%     % calculate CMatrix
%     eXVMatrix = [eXV_KRG eXV_SVR eXV_SHEP];
%     CMatrix   = srgtsWASComputeCMatrix(X, eXVMatrix)
% 
%     CMatrix =
% 
%     1.3144   -0.2050    0.2968
%    -0.2050    0.2229    0.0717
%     0.2968    0.0717    0.3766

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
% check inputs
if nargin == 2
    flag = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run
[NbPoints, NDV] = size(P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CMatrix
w = ones(NbPoints,1);

if flag == 1
    for c1 = 1 : NbPoints

        x = P(c1,:);

        for c2 = 1 : NDV
            if ( ( x(c2) == 0 ) || ( x(c2) == 1 ) )
                w(c1) = 0.5*w(c1);
            end
        end

    end
end

w       = diag(sparse(w));
CMatrix = ( (ErrorMatrix.')*w*ErrorMatrix )/NbPoints;

return
