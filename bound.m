function [ lu, opt ] = bound( func, n )
% get the lb and ub for benchmark function

switch func
    % TF1
    case 'sphere'         
        x_max=100;       x_min=-100;  opt = 0;
    case 'rosenbrock'         
        x_max=2.048;       x_min=-2.048;  opt = 0;
    case 'ackley'        
        x_max=32;       x_min=-32;  opt = 0;
    case 'griewank'         
        x_max=600;       x_min=-600;  opt = 0;
    case 'rastrigin'         
        x_max=5.12;       x_min=-5.12;  opt = 0;
        
    % TF
    case 'peak'
        x_min = -5;      x_max = 5;     opt = 0;
    case 'easom'
        x_min = 0;      x_max = 7;     opt = 0;
    case 'hart3'
        x_min = 0;      x_max = 1;     opt = 0;
    case 'shekel'
        x_min = 3;      x_max = 5;     opt = 0;
    case 'hart6'
        x_min = 0;      x_max = 1;     opt = 0;
    case 'ackley10'
        x_min = -0.6;      x_max = 0.6;     opt = 0;
        
    % cec benchmark
    case 'cec1'
        x_max=100;       x_min=-100;  opt = -450;
    case 'cec2'
        x_max=100;        x_min=-100;  opt = -450;
    case 'cec3'
        x_max=100;       x_min=-100;    opt = -450;
    case 'cec4'
        x_max=100;       x_min=-100;     opt = -450;
    case 'cec5'
        x_max=100;         x_min=-100;    opt = -310;
    case 'cec6'
        x_max=100;       x_min=-100;     opt = 390;
    case 'cec7'
        x_max=600;       x_min=0;     opt = -180;
    case 'cec8'
        x_max=32;        x_min=-32;     opt = -140;
    case 'cec9'
        x_max=5;         x_min=-5;       opt = -330;
    case 'cec10'
        x_max=5;         x_min=-5;       opt = -330;
    case 'cec11'
        x_max=0.5;      x_min=-0.5;     opt = 90;
    case 'cec12'
        x_max=pi;        x_min=-pi;      opt = -460;
    case 'cec13'
        x_max=1;         x_min=-3;       opt = -130;
    case 'cec14'
        x_max=100;       x_min=-100;     opt = 300;
    case 'cec15'
        x_max=2;       x_min=-5;     opt = 120;        
end 

lu = [x_min * ones(1, n); x_max * ones(1, n)];

end

