tic
options = optimoptions(@fminunc, 'Algorithm','quasi-newton', 'SpecifyObjectiveGradient', false);
options.Display = 'iter-detailed'; %'iter'
% options.SpecifyObjectiveGradient = true;
% options.Algorithm = 'quasi-newton'; %quasi-newton
options.CheckGradients = false;
% options.StepTolerance = 1e-6;
problem.options = options;

% problem.x0 = [+2.3400, -2.7400, +1.5800, +1.9500, +0.5000, -0.4000, -0.3000, +0.6000, +0.5500, +0.1400, +0.7400, +0.3700, +0.2800, +0.7500, +0.8500, +0.8900];
problem.x0 = [+1.1200, +1.2400, +2.4500, +2.1800, +0.5000, +0.4000, +0.7000, +0.5000, +0.3000, +0.8000, +0.6000, +0.4000, +0.5000, +0.7000, +0.7000, +0.3000];
% problem.x0 = [+2.4229, -2.4453, +1.8412,  +2.1660,  +0.2000, -0.0439, -0.3203, +0.5721, +0.5030, +0.2125, +0.8903, +0.4200, +0.1804, +0.7732, +0.8029, +0.6982];
problem.objective = @call_fx_m;
% problem.solver = 'fminunc';
[xo,fval,exitflag,output] = fminunc(problem);

fprintf('fx: %f. x: ', fval); fprintf('%8.4f ', xo); fprintf('\n');
toc