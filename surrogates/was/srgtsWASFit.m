function [srgtSRGT] = srgtsWASFit(srgtOPT)
%Function srgtsWASFit fits the specified weighted average surrogate.
% 
%    srgtSRGT = srgtsWASFit(srgtOPT)
% 
%srgtSRGT is the surrogate structure that contains the following fields:
%* NbSRGTs     : Number of basic surrogates.
%* SRGTS_ID    : Cell with the identifier of all basic surrogates.
%* srgtsSRGTs  : Cell with all basic surrogates.
%* WAS_Weights : Weights used in the linear combination of the surrogates.
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
%     srgtSRGTWAS = srgtsWASFit(srgtOPTWAS)
%
%     srgtSRGTWAS = 
% 
%         NbSRGTs: 3
%        SRGTS_ID: {'KRG'  'SVR'  'SHEP'}
%      srgtsSRGTs: {[1x1 struct]  [1x1 struct]  [1x1 struct]}
%     WAS_Weights: [0.0963 0.5677 0.3360]

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
srgtSRGT.NbSRGTs     = srgtOPT.NbSRGTs;
srgtSRGT.SRGTS_ID    = srgtOPT.SRGTS_ID;
srgtSRGT.srgtsSRGTs  = srgtOPT.srgtsSRGTs;
srgtSRGT.WAS_Weights = zeros(1, srgtOPT.NbSRGTs);

switch srgtOPT.WAS_Model
    case 'BestPRESS'
        PRESSRMS = srgtOPT.WAS_Options;
        [dummy, idx] = sort(PRESSRMS);
        idx = idx(1);
        srgtSRGT.WAS_Weights(idx) = 1;

    case {'OWSdiag', 'OWSconst', 'OWSfull'}
        CMatrix = srgtOPT.WAS_Options;
        if ( isequal(srgtOPT.WAS_Model, 'OWSdiag') || ...
                isequal(srgtOPT.WAS_Model, 'OWSconst') )
            Cinv = inv( diag( diag( CMatrix ) ) );
        else
            Cinv = pinv(CMatrix);
        end

        e      = ones(srgtOPT.NbSRGTs,1);
        lambda = 2/(e'*Cinv*e);

        srgtSRGT.WAS_Weights = (lambda*(Cinv*e)/2)';

        if isequal(srgtOPT.WAS_Model, 'OWSconst')
            A = [];
            b = [];
            Aeq = ones(size(srgtSRGT.WAS_Weights));
            beq = [1];
            lb = zeros(size(srgtSRGT.WAS_Weights.'));
            ub = [];

            srgtOPT = optimset('Display','off', ...
                'LargeScale', 'off');

            [x, jCWSoptc] = fmincon(@srgtsCWSObjectiveFunction,srgtSRGT.WAS_Weights.',A,b,Aeq,beq,lb,ub,[],srgtOPT,CMatrix);
            srgtSRGT.WAS_Weights = max(zeros(size(x)),x)';
        end


    case 'NPWS'
        PRESSRMS = srgtOPT.WAS_Options;
        for c1 = 1 : srgtOPT.NbSRGTs
            srgtSRGT.WAS_Weights(c1) = sum( PRESSRMS( [ [1:(c1 - 1)] [(c1 + 1) : end] ] ) );
        end

        srgtSRGT.WAS_Weights = srgtSRGT.WAS_Weights./( (srgtOPT.NbSRGTs - 1) * sum(PRESSRMS));

    case 'PWS'
        PRESSRMS = srgtOPT.WAS_Options{1};
        alpha    = srgtOPT.WAS_Options{2};
        beta     = srgtOPT.WAS_Options{3};

        Eavg = mean(PRESSRMS);
        wstar    = (PRESSRMS + alpha*Eavg).^beta;

        srgtSRGT.WAS_Weights = wstar./sum(wstar);

end
  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friend functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% objective function for MinCov
function f = srgtsCWSObjectiveFunction(x,C)
f = x'*C*x;
