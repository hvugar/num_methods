init_env

tic
options = optimoptions('fmincon');
options.Display = 'iter-detailed';
options.Algorithm = 'interior-point';
options.SpecifyObjectiveGradient = true;
options.CheckGradients = false;
options.StepTolerance       = 1.0000e-06;
options.OptimalityTolerance = 1.0000e-06;
options.FunctionTolerance   = 1.0000e-06;

problem.options = options;
problem.Aineq = []; 
problem.bineq = [];
problem.Aeq = [];
problem.beq = [];
problem.lb = [-1e+5, -1e+5, -1e+5, -1e+5, -1e+5, -1e+5, -1e+5, -1e+5, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]; 
problem.ub = [+1e+5, +1e+5, +1e+5, +1e+5, +1e+5, +1e+5, +1e+5, +1e+5, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95];

x0 = [];
[x0] = calllib('problem0H', 'init_strt_vector', x0);
%problem.x0 = x0;
%problem.x0 = [+1.1200, +1.2400, +1.4500, +1.1800, +0.5000, -0.4000, +0.7000, +0.5000, +0.4274, +0.6735, +0.6710, +0.3851, +0.5174, +0.7635, +0.5570, +0.4751];\
%problem.objective = @call_fx_m;
%problem.solver = 'fmincon';

%clc
%fprintf('Using algorithm: %s\n', options.Algorithm);
%fprintf('Initial point fx: %f. x0: ', call_fx_m(problem.x0)); fprintf('%8.4f ', problem.x0); fprintf('\n');

%[xo,fval,exitflag,output] = fmincon(problem);

%fprintf('fx: %f. x: ', fval); fprintf('%8.4f ', xo); fprintf('\n');
toc