function [recordBehaviour, flagRCD, sample] = VESAEA(benchmark, lu, sample, AdaptiveFE)
% This example demonstrates the powerful SPO algorithm
n = size(sample.x, 2);
fe = 0; flag = 1; 
flagRCD = [];
recordBehaviour = [];
while fe < AdaptiveFE
    
    f_min = min(sample.y);

    if flag       
        
        E_LOO = LOOCV(sample, f_min);
        srgtOPTKRG  = srgtsKRGSetOptions(E_LOO.x, E_LOO.y);
        LOO = srgtsKRGFit(srgtOPTKRG);        
        fhd = @(x)(-1*srgtsKRGEvaluate(x, LOO));      
        
        new_sample = PSO(fhd, n, lu);
        new_sample_y = benchmark(new_sample);
        fe = fe+1;
        sample.x = [sample.x; new_sample];
        sample.y = [sample.y; new_sample_y];
        
        if fe >= AdaptiveFE
            flagRCD = [flagRCD;flag];
            recordBehaviour = [recordBehaviour; min(sample.y) flag];
            break;            
        end
%         
        srgtOPTPRS  = srgtsRBFSetOptions(sample.x, sample.y);
        srgtSRGTPRS = srgtsRBFFit(srgtOPTPRS);
        fhd = @(x)(srgtsRBFEvaluate(x, srgtSRGTPRS));
       
%         [~, index] = sort(sample.y);
%         POP = sample.x(index, :);
%         POP = POP(1:2*n, :);
        new_sample = PSO(fhd, n, lu);
        new_sample_y = benchmark(new_sample);
        fe = fe+1;
        sample.x = [sample.x; new_sample];
        sample.y = [sample.y; new_sample_y];
        
    else

         
%         [lambda, gamma]=RBF(sample.x,sample.y, 'cubic');    
%         fhd=@(x) RBF_eval(x,sample.x,lambda,gamma,'cubic');
        
        % Monte Carlo simulation
        [~, index] = sort(sample.y);
        localSelect = max(floor(size(sample.x, 1)*0.1), 1);
%         localSelect = 1;
        sensitive = index(1:localSelect);
        
        
        srgtOPTRBF  = srgtsRBFSetOptions(sample.x, sample.y);
        srgtSRGTRBF = srgtsRBFFit(srgtOPTRBF);
        fhd = @(x)(srgtsRBFEvaluate(x, srgtSRGTRBF));
        
        new_sample = []; k = 0;
        while isempty(new_sample)
            new_sample = Monte_Carlo(sample, n, sensitive, lu, fhd);
            k = k+1;
            sensitive = index(1:localSelect+k);
        end
%         f_min = min(sample.y);
%         fhd = @(x)(criterion(srgtSRGTKRG, x, f_min, 'LCB', 0));
%         new_sample = PSO(fhd, n, lu);

        new_sample_y = benchmark(new_sample);  
        fe = fe+1;
        sample.x = [sample.x; new_sample];
        sample.y = [sample.y; new_sample_y];
    end
    flagRCD = [flagRCD;flag];
    recordBehaviour = [recordBehaviour; min(sample.y) flag];
    fprintf('VESAEA: sample index: %d Min: %f flag: %d \n', fe, min(sample.y), flag);
    if new_sample_y >= f_min
        flag = 1 - flag;  
    end

%     if mod(fe, n) == 0

%     end
end


sample.min = min(sample.y);
end

function sample = Monte_Carlo(sample, n, sensitive, lu, fhd)

NP = size(sample.x, 1) * n * 1000;
P_rand = repmat(lu(1, :), NP, 1) + rand(NP, n) .* (repmat(lu(2, :) - lu(1, :), NP, 1));
part = partition(P_rand, sample.x);
local_area = [];
for p=1:length(sensitive)
    local_area = [local_area; part{sensitive(p)}];
end

m = size(local_area, 1);
LOO_sensitive = zeros(m, 1);

for i=1:floor(m/10000)
    i1 = (i-1)*10000 + 1;
    i2 = i * 10000;
    LOO_sensitive(i1:i2) = feval(fhd, local_area(i1:i2, :));
end
r = mod(m, 10000);
LOO_sensitive(end-r+1:end) = feval(fhd, local_area(end-r+1:end, :));

[~, index_new] = min(LOO_sensitive);
sample = local_area(index_new, :);

end

function part = partition(P, c)

ds = pdist2(c, P);
[~, index] = min(ds, [], 1);
N = size(c, 1);
part = cell(1, N);
for i=1:N
    part{i} = P(find(index==i), :);
end

end


function f = local_search(sample, x, sensitive, LOO)

ds = pdist2(x, sample);
NP = size(x, 1);
f = zeros(NP, 1);
for i = 1:NP
    [~, index] = min(ds(i, :));
    if ismember(index, sensitive)
        f(i) = -1*srgtsKRGEvaluate(x(i, :), LOO);
    else
        f(i) = Inf;
    end
end

end


