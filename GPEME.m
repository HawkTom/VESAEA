function sample = GPEME(benchmark, lu, sample, AdaptiveFE)


n = size(sample.x, 2);

lambda = n; fe = 0;

while fe < AdaptiveFE
    %                 if mod(fe, 100)==0
%     fprintf('GPEME: sample index: %d Min: %f \n', fe, min(sample.y));
    %                 end
    [~, index] = sort(sample.y);
    pop = sample.x(index(1:lambda), :);
    fit = sample.y(index(1:lambda), :);
    
    offspring = zeros(lambda, n);
    for ind=1:lambda
        [~, index] = min(fit);
        gbest = pop(index, :);
        r = randi(lambda, 1, 2);
        v = gbest + 0.8 * (pop(r(1), :) - pop(r(2), :));
        jRand = randi(n);
        p1 = rand(1, n) <= 0.8;
        p1(jRand) = 1;
        u = pop(ind, :);
        u(p1) = v(p1);
        offspring(ind, :) = u;
    end
    
%     train.x = sample.x(end-tao+1:end, :);
%     train.y = sample.y(end-tao+1:end, :);
    train =sample;
    [train.x, ia, ~] = unique(train.x, 'rows');
    train.y = train.y(ia, :);
    srgtOPTPRS  = srgtsKRGSetOptions(train.x, train.y);
    KRG = srgtsKRGFit(srgtOPTPRS);
    [yhat1, mse1] = srgtsKRGPredictor(offspring, KRG);
    
    %                 [train.x, ia, ~] = unique(train.x, 'rows');
    %                 train.y = train.y(ia, :);
    %                 opts.type = 'BlindKriging';
    % %                 opts.regressionMaxOrder = 0;
    %                 srgtSRGTKRG = oodacefit( train.x, train.y, opts );
    %                 [yhat, mse] = srgtSRGTKRG.predict(offspring);
    appFit = yhat1 - 2 .* sqrt(mse1);
    
    %                 f_min = min(sample.y);
    %                 fhd = @(x)(criterion(KRG, x, f_min, 'LCB', 0));
    %                 appFit = feval(fhd, offspring);
    
    [~, index] = min(appFit);
    new_sample = offspring(index, :);
    new_sample_y = feval(benchmark, new_sample);
    fe = fe + 1;
    
    sample.x = [sample.x; new_sample];
    sample.y = [sample.y; new_sample_y];
end
sample.min = min(sample.y);
end









