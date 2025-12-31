% fdpendulum solver timings using timeit
% Example comparing odepsmh and odepsmJZ to ode45 for forced damped pendulum
% Richard Neidinger, 1/18/2021 to 2/28/25, 10/14/25 using timeit
format compact
format shortE

tspan = [0,200];
init = [0;2];  % y0 vector

stepsize = .6;  % for rk4 and odepsmh first setting
% ad in variable names stands for adaptive method
% ac in variable names stands for accurate adaptive method
% et in varoable names stands for elapsed time in sec

% make code for all psm solvers, find time of code generation
f4time = @() makepsmcode(@fdpendulum,2);
etcode = timeit(f4time)  % timeit averages across many runs
% this time is not included in timings for psm solvers

% rk4 Runge Kutta order 4
[trk4,yrk4] = rk4(@fdpendulum,tspan,stepsize,init);
f4time = @() rk4(@fdpendulum,tspan,stepsize,init); 
etrk = timeit(f4time,2)

% ode45 default options ('Refine' 4 interpolates four values for each step)
[tad,yad] = ode45(@fdpendulum,tspan,init);
f4time = @() ode45(@fdpendulum,tspan,init);
etad = timeit(f4time,2)

% ode45 with most accurate RelTol allowed
tolmin = 2.22045e-14;
options = odeset('RelTol',tolmin,'AbsTol',tolmin,'Refine',1);
[tac,yac] = ode45(@fdpendulum,tspan,init,options); 
f4time = @() ode45(@fdpendulum,tspan,init,options);
etac = timeit(f4time,2)

petac = 100; % percentage of ode45 accurate elapsed time
petrk = 100*etrk/etac;
petad = 100*etad/etac;
petcode = 100*etcode/etac;

% odepsmh
[tpsm,ypsm] = odepsmh(@fdpendulum,tspan,stepsize,init,20); 
f4time = @() odepsmh(@fdpendulum,tspan,stepsize,init,20); 
etpsm = timeit(f4time,2)
petpsm = 100*etpsm/etac;

% odepsmh with lower stepsize to match ode45 accurate result
stepsize2 = .25;
[tpsm2,ypsm2] = odepsmh(@fdpendulum,tspan,stepsize2,init,20); 
f4time = @() odepsmh(@fdpendulum,tspan,stepsize2,init,20); 
etpsm2 = timeit(f4time,2)
petpsm2 = 100*etpsm2/etac;

% odepsmJZ adaptive method with tolerance for high accuracy:
% tolmin = 2.22045e-14; same as for smallest RelTol allowed by ode45
[tpsmac,ypsmac,degac] = odepsmJZ(@fdpendulum,tspan,init,tolmin);
f4time = @() odepsmJZ(@fdpendulum,tspan,init,tolmin);
etpsmac = timeit(f4time,3)
petpsmac = 100*etpsmac/etac;

% odepsmJZ adaptive method with tolerance for graphic viewing accuracy
tol = 1e-3;
[tpsmad,ypsmad,degad] = odepsmJZ(@fdpendulum,tspan,init,tol);
f4time = @() odepsmJZ(@fdpendulum,tspan,init,tol);
etpsmad = timeit(f4time,3)
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
fprintf('One-time generation of series code for psm solvers takes rel time %4.1f%%.\n',petcode)