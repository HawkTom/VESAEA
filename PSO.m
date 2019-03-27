function gbest = PSO(fhd, c, lu)
% Usage: [ time,gbest,POP ] =PSO_D(Data,bu,bd )
% -----------------------------------------------------------------
% Important Note: This code needs intralled SURROGATE TOOLBOX(https://sites.google.com/site/srgtstoolbox/)
% -----------------------------------------------------------------
% Input:
% Data          - Data with c Decision Variables and Exact Objective Value
% bu            - Upper Boundary of c Decision Variables
% bd            - Lower Boundary of c Decision Variables
%
% Output: 
% time          - Execution Time
% Record        - Logbook of Evaluated Solutions
% gbest         - Predicted Optimum over Generations with c Decision Variables
%
    %%%%    Authors:    Handing Wang, Yaochu Jin, John Doherty
    %%%%    University of Surrey, UK
    %%%%    EMAIL:      wanghanding.patch@gmail.com
    %%%%    WEBSITE:    https://sites.google.com/site/handingwanghomepage
    %%%%    DATE:       May 2018
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:
%Handing Wang, Yaochu Jin, John Doherty, Committee-based Active Learning for Surrogate-Assisted Particle Swarm Optimization of Expensive Problems, IEEE Transactions on Cybernetics, vol.47, no.9, pp.2664-2677, 2017.
%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------
bu = lu(2, 1);
bd = lu(1, 1);

n=50;
POP = initialize_pop(n,c,bu,bd);
obj = feval(fhd, POP);

POP=[POP,obj];
lbest=POP;

[best,Ib]=min(POP(:,c+1));
gbest=POP(Ib,:);
v=(2*rand(n,c)-1).*(ones(n,1)*(bu-bd))*0.5;
g=0;
gmax=200;
B=gbest;
while g<gmax  
    [ POP,v ] = fly(POP,bu,bd,gbest,lbest,v,g/gmax);
    obj = feval(fhd, POP(:, 1:c));

    POP(:,c+1)=obj;
    [ POP,gbest,lbest] = Update_best(POP,gbest,lbest);
    best=gbest(end);
    if best<B(end,end)
        B=[B;gbest];
    end
    g=g+1;
end

gbest = gbest(:, 1:c);

end
