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

% problem.x0 = [+2.3400, -2.7400, +1.5800, +1.9500, +0.5000, -0.4000, -0.3000, +0.6000, +0.5500, +0.1400, +0.7400, +0.3700, +0.2800, +0.7500, +0.8500, +0.8900];
problem.x0 = [+1.1200, +1.2400, +1.4500, +1.1800, +0.5000, -0.4000, +0.7000, +0.5000, +0.4274, +0.6735, +0.6710, +0.3851, +0.5174, +0.7635, +0.5570, +0.4751];
% problem.x0 = [+2.4229, -2.4453, +1.8412,  +2.1660,  +0.2000, -0.0439, -0.3203, +0.5721, +0.5030, +0.2125, +0.8903, +0.4200, +0.1804, +0.7732, +0.8029, +0.6982];
% problem.x0 = [+0.5636, -0.2421, +0.2505, -0.7679, +0.3220, -0.6179, -0.1833, -0.3160, +0.4269, +0.3364, +0.8651, 0.2336, 0.5436, 0.5788, 0.5269, 0.5860];

problem.objective = @call_fx_m;
problem.solver = 'fmincon';

clc
fprintf('Using algorithm: %s\n', options.Algorithm);
fprintf('Initial point fx: %f. x0: ', call_fx_m(problem.x0)); fprintf('%8.4f ', problem.x0); fprintf('\n');

[xo,fval,exitflag,output] = fmincon(problem);

fprintf('fx: %f. x: ', fval); fprintf('%8.4f ', xo); fprintf('\n');
toc