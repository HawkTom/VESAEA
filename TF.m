function y = TF1(X, name)
% test function for optimization
[pop_size, fun_dim] = size(X);
switch name
    % Sphere
    case 'sphere'
        bias = [-40,-13,-15,-40,-54,-70,56,31,-36,-16,-31,-10,-36,17,46,-75,-33,76,32,-47,-8,58,33,31,45,46,48,49,-20,7];
        X = X - bias(1:fun_dim);
        y = sum( (X) .^ 2, 2);
    % Rosenbrock
    case 'rosenbrock'
        bias = [-2,0,0,-1,-1,1,-2,-1,-1,1,-1,0,-1,0,-1,0,-1,-1,0,-2,-1,0,-1,0,1,-1,1,-2,-1,-1];
        X = X - bias(1:fun_dim);
        y = 100.0 * sum(((X(:, 1 : (fun_dim - 1)) .^ 2) - X(:, 2 : fun_dim)) .^ 2, 2) + ...
            sum((X(:, 1 : (fun_dim - 1)) - 1) .^ 2, 2);
    % Ackley
    case 'ackley'
        bias = [-4,22,16,19,-1,-15,23,-18,-2,-22,-3,17,-3,20,-4,-6,9,-9,10,-17,-4,6,-1,16,-4,-10,-9,2,-23,-18];
        X = X - bias(1:fun_dim);
        a=20;
        b=0.2;
        cc=2*pi;
        y=0-a.*exp(0-b.*(mean((X).^2,2)).^0.5)-exp(mean(cos(cc.*(X)),2))+a+exp(1);
    % Griewank
    case 'griewank'        
        bias = [-189,296,-313,-192,4,-13,104,134,-311,479,16,-313,2,-172,147,313,-394,-305,-9,-420,167,-38,-115,71,80,327,-349,323,-85,385];
        X = X - bias(1:fun_dim);
        y = sum(X .^ 2, 2) / 4000.0 - ...
            prod(cos((X) ./ sqrt(repmat(1 : fun_dim, pop_size, 1))), 2) + 1;
    % Rastrigin
    case 'rastrigin'
        bias = [-2,2,3,2,-3,-1,0,2,-3,-1,-4,-1,1,-5,0,2,1,-1,-2,1,0,-2,3,-2,2,-4,2,-1,2,-4];
        X = X - bias(1:fun_dim);
        y = sum((X) .^ 2 - 10.0 * cos(2.0 * pi * (X)) + 10, 2);
end

end