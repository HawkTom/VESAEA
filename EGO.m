function sample = EGO(benchmark, lu, sample, AdaptiveFE)

    n = size(sample.x, 2);
    %% Main Loop
    for i=1:AdaptiveFE
        % kriging
        srgtOPTKRG  = srgtsKRGSetOptions(sample.x, sample.y);
        srgtSRGTKRG = srgtsKRGFit(srgtOPTKRG);
        f_min = min(sample.y);
        fhd = @(x) criterion(srgtSRGTKRG, x, f_min, 'LCB', 0);
        new_sample = cmaes(fhd, n, lu);
        sample.x = [sample.x; new_sample];
        sample.y = [sample.y; benchmark(new_sample)];
%         plot3(new_sample(1, 1), new_sample(1, 2), benchmark(new_sample), '.', 'Color', 'r', 'MarkerSize', 20)
%         fprintf('EGO: sample index: %d Min: %f\n', i, min(sample.y))
    end
    sample.min = min(sample.y);
end


function obj = criterion(KRG, x, f_min, option, iter)
%calculate expected improvment
[y, mse] = srgtsKRGPredictor(x,KRG);
s = sqrt(max(0, mse));

switch option
    
    case 'ExI'
        % Expected Improvement
        ei = (f_min - y).*normcdf((f_min - y)./s) + s.*normpdf((f_min - y)./s);
        obj = -ei;
    case 'LCB'
        % Lower Confidence Bound
        lcb = y - 2 * s;
        obj = lcb;
    case 'Mean'
        % Mean
        obj = y;
    case 'MSE'
        % uncertain
        obj = s;
    case 'PoI'
        % Probability of Improvement
        poi = normcdf((y - f_min)./s);
        obj = poi;
    case 'MIX'
        p = gridsamp([0.1; 0.9], 400);
        if rand() > p(iter)
            obj = normcdf((f_min - y)./s);
        else
            ei = (f_min - y).*normcdf((f_min - y)./s) + s.*normpdf((f_min - y)./s);
            obj = -ei;
        end
    otherwise
        error('No such option')
end
end
