function predvar = srgtsWASPredictionVariance(x, srgtSRGT)
%Function srgtsWASPredictionVariance returns the estimated prediction
%variance of a weighted average surrogate. Thus, for example:
%
%     PREDVAR = srgtsWASPredictionVariance(X, SURROGATE): returns the
%     prediction variance PREDVAR estimated by the weighted average
%     surrogate SURROGATE at all X sites. X can be either a single row
%     vector (single point, with each column representing a variable) or a
%     matrix (each row represents a point). PREDVAR is an NPOINTS-by-1
%     vector.
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
%     % weighted average surrogate
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
%     CMatrix   = srgtsWASComputeCMatrix(X, eXVMatrix);
%
%     % options for WAS
%     srgtsOPTs   = {srgtOPTKRG  srgtOPTSVR  srgtOPTSHEP};
%     srgtsSRGTs  = {srgtSRGTKRG srgtSRGTSVR srgtSRGTSHEP};
%     WAS_Model   = 'OWSdiag';
%     WAS_Options = CMatrix;
% 
%     srgtOPTWAS  = srgtsWASSetOptions(srgtsOPTs, srgtsSRGTs, WAS_Model, WAS_Options);
%     srgtSRGTWAS = srgtsWASFit(srgtOPTWAS);
% 
%     % create test points
%     Xtest = linspace(designspace(1), designspace(2), 100)';
%
%     % evaluate surrogate at Xtest
%     PredVar = srgtsWASPredictionVariance(Xtest, srgtSRGTWAS);
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

npoints = length(x(:,1));
predvar = zeros(npoints,1);
for c1 = 1 : npoints

    ysrgts = zeros(srgtSRGT.NbSRGTs, 1);
    for c2 = 1 : srgtSRGT.NbSRGTs
        eval(sprintf('ysrgts(c2) = srgts%sEvaluate(x(c1,:), srgtSRGT.srgtsSRGTs{c2});', srgtSRGT.SRGTS_ID{c2}));
    end
    
    predvar(c1) = std(ysrgts)^2;
    
end

return
