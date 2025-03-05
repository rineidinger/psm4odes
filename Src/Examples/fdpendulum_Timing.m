% fdpendulum solver timings repeating for numreps
% Example comparing odepsmh and odepsmJZ to ode45 for forced damped pendulum
% Richard Neidinger, 1/18/2021 to 2/28/25
numreps = 10;  % will average runtime over numreps calls of each solver
format compact
format shortE

tspan = [0,200];
init = [0;2];  % y0 vector

stepsize = .6;  % for rk4 and odepsmh first setting
% ad in variable names stands for adaptive method
% ac in variable names stands for accurate adaptive method
% et in varoable names stands for elapsed time in sec

% rk4 Runge Kutta order 4
tic; for i=1:numreps
    [trk4,yrk4] = rk4(@fdpendulum,tspan,stepsize,init); 
end
etrk = toc/numreps 

% ode45 default options ('Refine' 4 interpolates four values for each step)
tic; for i=1:numreps
    [tad,yad] = ode45(@fdpendulum,tspan,init); 
end
etad = toc/numreps

% ode45 with most accurate RelTol allowed
tolmin = 2.22045e-14;
options = odeset('RelTol',tolmin,'AbsTol',tolmin,'Refine',1);
tic; for i=1:numreps
    [tac,yac] = ode45(@fdpendulum,tspan,init,options); 
end
etac = toc/numreps
petac = 100; % percentage of ode45 accurate elapsed time
petrk = 100*etrk/etac;
petad = 100*etad/etac;

% odepsmh
makepsmcode(@fdpendulum,2); % do not include one-time writing series code
tic; for i=1:numreps
    [tpsm,ypsm] = odepsmh(@fdpendulum,tspan,stepsize,init,20); 
end
etpsm = toc/numreps
petpsm = 100*etpsm/etac;

% odepsmh with lower stepsize to match ode45 accurate result
stepsize2 = .25;
tic; for i=1:numreps
    [tpsm2,ypsm2] = odepsmh(@fdpendulum,tspan,stepsize2,init,20); 
end
etpsm2 = toc/numreps
petpsm2 = 100*etpsm2/etac;

% odepsmJZ adaptive method with tolerance for high accuracy:
% tolmin = 2.22045e-14; same as for smallest RelTol allowed by ode45
tic; for i=1:numreps
    [tpsmac,ypsmac,degac] = odepsmJZ(@fdpendulum,tspan,init,tolmin);
end
etpsmac = toc/numreps
petpsmac = 100*etpsmac/etac;

% odepsmJZ adaptive method with tolerance for graphic viewing accuracy
tol = 1e-3;
tic; for i=1:numreps
    [tpsmad,ypsmad,degad] = odepsmJZ(@fdpendulum,tspan,init,tol);
end
etpsmad = toc/numreps
petpsmad = 100*etpsmad/etac;

% adaptive stepsizes
dtac = tac(2:end) - tac(1:end-1);  % dt = h's for accurate ode45
dtad = tad(2:end) - tad(1:end-1);  % dt of points for default ode45
dtad = 4*dtad(1:4:end); % actual since ode45 uses Refine 4 by default (see Help on odeset)
dtpsmad = tpsmad(2:end) - tpsmad(1:end-1);
dtpsmac = tpsmac(2:end) - tpsmac(1:end-1);

% Find true (best) ending value using Symbolic Toolbox
% output from run below to avoid running this section of code
% and tspan = [0,200]; init = [0;2];
best = 17.417045282416465580661324578379; % with stepsize .25
%      17.417045282416465582192310843963 with stepsize .125
% using vpa (Variable Precision Arithmetic default 32 digits) 
% using odepsmh to degree 25 but it takes about 8 minutes on my old laptop
% stepsize2 = .25;
% tic; [tpsmvh,ypsmvh] = odepsmh(@fdpendulum,tspan,vpa(stepsize2),vpa(init),25); etpsmvh = toc;
% best = ypsmvh(end,1)  % assumed to be the most accurate value
% note: automatic conversion vpa(.1) from DE keeps the quadruple precision

% % try odepsmJZ with vpa; faster and reliable
% tol2 = 1e-20;
% tic; [tpsmv,ypsmv,degv] = odepsmJZ(@fdpendulum,tspan,vpa(init),tol); etpsmv = toc;
% dtpsmv = tpsmv(2:end) - tpsmv(1:end-1);
% best = ypsmv(end,1);  % assumed to be the most accurate value

% Table of results:
fprintf('%23s %10s %20s %8s %9s\n','num steps','mean step','y(end) value','error','rel time');
str = ': %6d %10.4f %20.14f %8.1e %8.1f%%\n';  % time in sec could use 8.1e or 8.5f
fprintf(['  runge kutta 4',str],length(trk4)-1, stepsize, yrk4(end,1),abs(yrk4(end,1)-best),petrk);
fprintf(['  ode45 default',str],length(dtad), mean(dtad), yad(end,1),abs(yad(end,1)-best),petad);
fprintf([' ode45 accurate',str],length(dtac), mean(dtac), yac(end,1),abs(yac(end,1)-best),petac);
fprintf([' odepsmh deg 20',str],length(tpsm)-1, stepsize, ypsm(end,1),abs(ypsm(end,1)-best),petpsm);
fprintf([' odepsmh deg 20',str],length(tpsm2)-1, stepsize2, ypsm2(end,1),abs(ypsm2(end,1)-best),petpsm2);
fprintf(['odepsmJZ deg %2d',str],degac,length(dtpsmac), mean(dtpsmac), ypsmac(end,1),abs(ypsmac(end,1)-best),petpsmac);
fprintf(['odepsmJZ deg %2d',str],degad,length(dtpsmad), mean(dtpsmad), ypsmad(end,1),abs(ypsmad(end,1)-best),petpsmad);