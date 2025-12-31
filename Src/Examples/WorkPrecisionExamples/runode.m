function [yend,scd,time,numsteps,deg] = runode(solverhandle,dehandle,tspan,init,options,refsoln,once)
% Testing time and accuracy (given refsoln) of ode solver.
% once true uses toc-toc; once false uses timeit to average over multiple runs.
% Richard Neidinger 12/6/25
solvername = func2str(solverhandle);
dename = func2str(dehandle);
deg = NaN;  % only computed for odepsmJZ
if nargin < 7, once = true; end

if strcmp(solvername,'odepsmJZ')
    options = odeget(options,'RelTol');
    tic;
    [t,y,deg] = feval(solverhandle,dehandle,tspan,init,options);
    time = toc;
    solvername = [solvername,' ',num2str(deg)];
else
    tic;
    [t,y] = feval(solverhandle,dehandle,tspan,init,options);
    time = toc;
end
dt = t(2:end) - t(1:end-1);
numsteps = length(t)-1;
compnum = length(refsoln);
yend = y(end,1:compnum);
scd = -log10(norm((yend-refsoln)./refsoln,Inf));

if ~once  % not once
    f4time = @() feval(solverhandle,dehandle,tspan,init,options);
    time = timeit(f4time,2);
    % anonymous functions have only 1 output; second arg of timeit for multiple
end

% Table of results:
fprintf('%-12s: %13s: %9s %10s %10s %20s %8s %10s\n','solver','de','steps','min step','mean step','y1(end) value','scd','seconds');
str = '%12s: %13s: %9d %10.1e %10.1e %22.15e %6.1f %10.4f\n'; 
fprintf(str,solvername,dename,numsteps,min(dt),mean(dt),yend(1),scd,time);