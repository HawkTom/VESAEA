clear; close all
warning('off');
addpath(genpath('surrogates'))
addpath('CALSAPSO')

% rng(20181010);
% Problem Definition
Dim = [5, 10, 15, 20, 30];
name = {'sphere', 'rosenbrock', 'ackley', 'griewank', 'rastrigin'};
start = tic;
for run=1:25
    for i = 1:5 % iterate for dimension
        for j =1:5  % iterate for function
            
            n = Dim(i);
            func_name = name{j};
            lu = bound(func_name, n);
            benchmark = @(x)(TF(x, func_name));
            
            MAX_FE = 10;
            InitialFE = 2*n;
            AdaptiveFE =3*n;
            
%             Samples for initilization with LHD
            sample.x = repmat(lu(1, :), InitialFE, 1) + lhsdesign(InitialFE, n, 'iterations', 1000) .* (repmat(lu(2, :) - lu(1, :), InitialFE, 1));
            sample.y = benchmark(sample.x);
            
            
            %             Funciton Interface
            %             Input:
            %                 banchmark: the handle of test function
            %                 lu: the search bound of test function
            %                 sample: initial samples
            %                 AdaptiveFE: number of function evalution
            %             Ouput:
            %                 sample: the final sample of algorithms
            result = cell(1, 5);
            result{1} = VESAEA(benchmark, lu, sample, AdaptiveFE);            
            result{2} = EGO(benchmark, lu, sample, AdaptiveFE);
            result{3} = GPEME(benchmark, lu, sample, AdaptiveFE);
            result{4} = SSLAPSO(benchmark, lu, sample, AdaptiveFE);

            Data = [sample.x sample.y];
            BU = lu(2, :);
            BD = lu(1, :);
            new_samples = CALSAPSO( Data,BU,BD,benchmark, AdaptiveFE);
            result{5}.x = [sample.x; new_samples(:, 1:n)];
            result{5}.y = [sample.y; new_samples(:, n+1)];
            result{5}.min = min(result{5}.y);
            
            
            filename = strcat('result/result_run',num2str(run),'_', name{j}, '_', num2str(n), '.mat');
%             save(filename, 'result', 'func_name')
            fprintf(strcat('\n', filename, '\n'));
            str = sprintf('RUN:%d, FUN:%d, DIM: %d',run,j, n);
            fprintf(str);
         end
    end
end
toc(start)









