function [ E_LOO ] = LOOCV( sample, f_min )
% valid-one-out cross validation

N = size(sample.x, 1);
E_LOO.x = sample.x;
E_LOO.y = zeros(N, 1);

for i=1:N
    valid.x = sample.x(i, :);
    valid.y = sample.y(i);
    
    train =sample;
    train.x(i, :) = [];
    train.y(i) = [];  
    
    % kriging
%     srgtOPTKRG  = srgtsKRGSetOptions(train.x, train.y);
%     srgtSRGTKRG = srgtsKRGFit(srgtOPTKRG);
%     
%     yhat = srgtsKRGEvaluate(valid.x, srgtSRGTKRG);

    srgtOPTPRS  = srgtsPRSSetOptions(train.x, train.y);
    srgtSRGTPRS = srgtsPRSFit(srgtOPTPRS);
    yhat = srgtsPRSEvaluate(valid.x, srgtSRGTPRS);
    
%     opts.type = 'BlindKriging';
%     opts.regressionMaxOrder = 0;
%     srgtSRGTKRG = oodacefit( train.x, train.y, opts );
%     yhat = srgtSRGTKRG.predict(valid.x);
%     E_LOO.y(i) = abs(valid.y - yhat) / (valid.y - 0.5*f_min);
    E_LOO.y(i) = abs(valid.y - yhat);
    
end

end

