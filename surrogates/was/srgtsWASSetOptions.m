function srgtOPT = srgtsWASSetOptions(srgtsOPTs, srgtsSRGTs, WAS_Model, WAS_Options)
%Function srgtsWASSetOptions creates the SURROGATES Toolbox option
%structure for weighted average surrogate. This structure contains the
%following fields:
%
%* GENERAL PARAMETERS
%
%   SRGT        - Identifier of the surrogate technique: 'WAS'.
%
%* WEIGHTED AVERAGE SURROGATE PARAMETERS
%
%   WAS_Model   - Strategy used for weighting selection [ string |
%                 'BestPRESS' | 'OWSdiag' | 'OWSconst' | 'OWSfull' | 'PWS' |
%                 'NPWS' ]. Default: 'BestPRESS'.
%                   * BestPRESS : surrogate with smallest PRESS value.
%                   * OWSdiag   : Optimal weighted surrogates (minimization
%                                 of mean square error), but using only the
%                                 diagonal elements of the C-matrix (matrix
%                                 of the mean square errors).
%                   * OWSconst  : Optimal weighted surrogates, using the
%                                 full C-matrix but with a positivity
%                                 constraint in the optimization.
%                   * OWSfull   : Optimal weighted surrogates, using the
%                                 full C-matrix.
%                   * PWS       : PRESS weighted surrogate.
%                   * NPWS      : Non-parametric PRESS weighted surrogate.
%
%   WAS_Options - WAS_Model paramenters. Therefore, there is no default
%                 value. Here it is the description for each WAS_Model:
%                   * BestPRESS: 1 x NbSurrogates vector, which elements
%                                are the root mean square PRESS
%                                (PRESSrms) of each surrogate.
%                   * OWSdiag  : C-matrix, i.e. the matrix of mean
%                                square errors.
%                   * OWSconst : C-matrix.
%                   * OWSfull  : C-matrix.
%                   * PWS      : {p1 p2 p3},  where p1, p2 and p3 are
%                                elements of a cell array. Parameter p1 is
%                                the 1 x NbSurrogates vector, which
%                                elements are the root mean square PRESS
%                                (PRESSrms) of each surrogate. Parameters
%                                p2 and p3 are the alpha and beta
%                                parameters, which control the importance
%                                of the averaging and the importance of
%                                individual surrogate, respectively (as a
%                                suggestion, p2 = 0.05 and p3 = -1).
%                   * NPWS     : 1-by-NBSURROGATES vector, which elements
%                                are the root mean square PRESS (PRESSrms)
%                                of each surrogate.
%
%   NbSRGTs     - Number of basic surrogates.
%
%   srgtsSRGTs  - Cell with all basic surrogates.
%
%   SRGTS_ID    - Cell with the identifier of all basic surrogates.
%
%This is how you can use srgtsWASSetOptions:
% 
%     OPTIONS = srgtsWASSetOptions: creates a structure with the empty
%     parameters.
%
%     OPTIONS = srgtsWASSetOptions(P, T): Given the sampled data P (input
%     variables) and T (output variables), it creates a structure with
%     default parameters used for all not specified fields.
%
%     OPTIONS = srgtsWASSetOptions(P, T, ...
%         WAS_Kernel, WAS_KernelOptions, WAS_C, WAS_Loss,
%         WAS_Insensitivity)): it creates a structure with each of the
%         specified fields.
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
%     srgtOPTWAS  = srgtsWASSetOptions(srgtsOPTs, srgtsSRGTs, WAS_Model, WAS_Options)
% 
%     srgtOPTWAS = 
% 
%            SRGT: 'WAS'
%       WAS_Model: 'OWSdiag'
%     WAS_Options: [3x3 double]
%         NbSRGTs: 3
%      srgtsSRGTs: {[1x1 struct]  [1x1 struct]  [1x1 struct]}
%        SRGTS_ID: {'KRG'  'SVR'  'SHEP'}

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
% options
srgtOPT.SRGT = 'WAS';
srgtOPT.WAS_Model   = WAS_Model;
srgtOPT.WAS_Options = WAS_Options;

srgtOPT.NbSRGTs    = length(srgtsSRGTs);
srgtOPT.srgtsSRGTs = srgtsSRGTs;
for c1 = 1 : srgtOPT.NbSRGTs
    srgtOPT.SRGTS_ID{c1} = srgtsOPTs{c1}.SRGT;
end

return
