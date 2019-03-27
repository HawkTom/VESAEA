function DB = SSLAPSO(benchmark, lu, sample, AdaptiveFE)

n = size(sample.x, 2);

pop_size = min(50, 2*n);
[~, index] = sort(sample.y);
pop = sample.x(index(1:pop_size), :);
fit = sample.y(index(1:pop_size), :);
v = zeros(pop_size, n);

pbest = pop; pbestf = fit;

gbestf = fit(1);
gbest = pop(1, :);

fe = 0;
DB = sample;
EDB = DB;
while fe < AdaptiveFE
    xBar = mean(pop);
    [pbestf, index] = sort(pbestf);
    pbest = pbest(index, :);
    pop = pop(index, :);
    fit = fit(index);
    for ind=1:pop_size
        r = rand(1, 3);
        ps = randi([ind, pop_size]);
        v(ind, :) = r(1) * v(ind, :) + r(2) * (pbest(ps, :) - pop(ind, :)) + r(3) * 0.001 * (xBar - pop(ind, :));
    end
    pop = pop + v;
    %                 for ind = 1:pop_size
    %                     pop(ind,:) = max(pop(ind,:), lu(1,:));
    %                     pop(ind,:) = min(pop(ind,:), lu(2,:));
    %                 end
    %                 train_DB.x = DB.x(end-2*n-1:end, :);
    %                 train_DB.y = DB.y(end-2*n-1:end, :);
    train_DB = DB;
    [train_DB.x, ia, ~] = unique(train_DB.x, 'rows');
    train_DB.y = train_DB.y(ia, :);
    OPT  = srgtsRBFSetOptions(train_DB.x, train_DB.y);
    RBF_DB = srgtsRBFFit(OPT);
    
    
    %                 train_EDB.x = EDB.x(end-4*n-3:end, :);
    %                 train_EDB.y = EDB.y(end-4*n-3:end, :);
    train_EDB = EDB;
    [train_EDB.x, ia, ~] = unique(train_EDB.x, 'rows');
    train_EDB.y = train_EDB.y(ia, :);
    OPT  = srgtsRBFSetOptions(train_EDB.x, train_EDB.y);
    RBF_EDB = srgtsRBFFit(OPT);
    
    
    appFit_DB = srgtsRBFEvaluate(pop, RBF_DB);
    appFit_EDB = srgtsRBFEvaluate(pop, RBF_EDB);
    
    % update the personal best position
    flag = 1;
    for ind=1:pop_size
        if appFit_DB(ind) < pbestf(ind) && appFit_EDB(ind) < pbestf(ind)
            fit(ind) = feval(benchmark, pop(ind, :));
            DB.x = [DB.x; pop(ind, :)];
            DB.y = [DB.y; fit(ind)];
            EDB.x = [EDB.x; pop(ind, :)];
            EDB.y = [EDB.y; fit(ind)];
            fe = fe + 1;
            flag = 0;
        else
            fit(ind) = max(appFit_DB(ind), appFit_EDB(ind));
        end
        if fit(ind) < pbestf(ind)
            pbestf(ind) = fit(ind);
            pbest(ind, :) = pop(ind, :);
        end
    end
    
    % select unlabelled data
    appIndex = find(feval(benchmark, pop)~= fit);
    appPOP = pop(appIndex, :);
    appFIT = fit(appIndex);
    
    trueIndex = find(feval(benchmark, pop) == fit);
    truePOP = pop(trueIndex, :);
    trueFIT = fit(trueIndex);
    
%     compare = [feval(benchmark, pop) fit];
    if ~isempty(trueFIT)
        UL = selectUnlabeled(EDB, appPOP, appFIT, truePOP, trueFIT, n);
        EDB.x = [EDB.x; UL.x];
        EDB.y = [EDB.y; UL.y];
    end
    
    
    % update the global best position
    [~, k] = min(pbestf);
    
    if pbestf(k) < gbestf
        if pbestf(k) ~= feval(benchmark, pbest(k, :))
            pbestf(k) = feval(benchmark, pbest(k, :));
            fe = fe + 1;
            DB.x = [DB.x; pbest(k, :)];
            DB.y = [DB.y; pbestf(k)];
            EDB.x = [EDB.x; pbest(k, :)];
            EDB.y = [EDB.y; pbestf(k)];
        end
        if pbestf(k) < gbestf
            gbestf = pbestf(k);
            gbest = pbest(ind, :);
        end
    else
        if flag
            [~, min_index] = min(fit);
            fit(min_index) = feval(benchmark, pop(min_index, :));
            fe = fe + 1;
            if fit(min_index) < pbestf(min_index)
                pbestf(min_index) = fit(min_index);
                pbest(min_index, :) = pop(min_index, :);
            end
            DB.x = [DB.x; pop(min_index, :)];
            DB.y = [DB.y; fit(min_index)];
            EDB.x = [EDB.x; pop(min_index, :)];
            EDB.y = [EDB.y; fit(min_index)];
            if fit(min_index) < gbest
                gbest = pop(min_index, :);
                gbestf = fit(min_index);
            end
        end
    end
    
%     fprintf('GPEME: sample index: %d Min: %f \n', fe, min(DB.y));
    
end
DB.min = min(DB.y);
end


function UL = selectUnlabeled(EDB, pop1, fit1, pop2, fit2, n)

%     train.x = EDB.x(end-4*n-3:end, :);
%     train.y = EDB.y(end-4*n-3:end, :);
train = EDB;
[train.x, ia, ~] = unique(train.x, 'rows');
train.y = train.y(ia, :);
OPT  = srgtsRBFSetOptions(train.x, train.y);
RBF_EDB = srgtsRBFFit(OPT);
appEDB = srgtsRBFEvaluate(pop2, RBF_EDB);

FAI = zeros(size(fit1));
for i = 1:length(fit1)
    tmp.x = [train.x; pop1(i, :)];
    tmp.y = [train.y; fit1(i, :)];
    [tmp.x, ia, ~] = unique(tmp.x, 'rows');
    tmp.y = tmp.y(ia, :);
    OPT  = srgtsRBFSetOptions(tmp.x, tmp.y);
    RBF_TMP = srgtsRBFFit(OPT);
    appTMP = srgtsRBFEvaluate(pop2, RBF_TMP);
    
    FAI(i) = min(abs(appTMP - fit2) - abs(appEDB - fit2));
end
[~, index] = max(FAI);
UL.x = pop1(index, :);
UL.y = fit1(index, :);
end






