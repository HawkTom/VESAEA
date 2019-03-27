clear; close all
warning('off');
addpath(genpath('surrogates'))
addpath('CALSAPSO')

% random seed
% rng(20190327);


% Problem Definition
Dim = [5, 10, 15, 20, 30];
% name = {'peak', 'easom', 'hart3', 'shekel', 'hart6', 'ackley10'};
name = {'sphere', 'rosenbrock', 'ackley', 'griewank', 'rastrigin'};
start = tic;
for run=1:25
    for i = 1:5 
        for j = 1:5 % iterate for function
            
            
            n = Dim(i);
            func_name = name{j};
            lu = bound(func_name, n);
            benchmark = @(x)(TF(x, func_name));
            
            
            InitialFE = 2*n;
            AdaptiveFE = 3*n;
            
            
%          Samples for initilization with LHD
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
            
            result = cell(1, 2);

            [recordMin, flagRCD, result{1}] = VESAEA(benchmark, lu, sample, AdaptiveFE);  
            result{1}.flag = flagRCD;
            result{1}.recordMin = recordMin;
            
            
            [recordMin, flagRCD, result{2}] =VESAEA_woVLS(benchmark, lu, sample, AdaptiveFE);  
            result{2}.flag = flagRCD;
            result{2}.recordMin = recordMin;

                
            filename = strcat('result/gl_result_run',num2str(run),'_', name{j}, '_', num2str(n), '.mat');
%             save(filename, 'result', 'func_name')
            fprintf(strcat('\n', filename, '\n'));
            str = sprintf('RUN:%d, compare: 1. %f,  2. %f\n',run, result{1}.min, result{2}.min);
            fprintf(str);
         end
    end
end
toc(start)









