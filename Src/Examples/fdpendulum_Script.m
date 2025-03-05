% Example comparing odepsmh and odepsmJZ to ode45 for forced damped pendulum
% Richard Neidinger, 1/18/2021 to 3/5/25
format compact
format shortE

tspan = [0,200];
init = [0;2];  % y0 vector

stepsize = .6;  % for rk4 and odepsmh
% ad in variable names stands for adaptive method
% ac in variable names stands for accurate adaptive method
% et in varoable names stands for elapsed time in sec

% Use fdpendulum_Timing to see relative timings 
% The one-run timings of this script are unreliable and not reported

% rk4 Runge Kutta order 4
tic; [trk4,yrk4] = rk4(@fdpendulum,tspan,stepsize,init); etrk = toc;

% ode45 default options ('Refine' 4 interpolates four values for each step)
tic; [tad,yad] = ode45(@fdpendulum,tspan,init); etad = toc;

% ode45 with most accurate RelTol allowed by ode45
tol = 2.22045e-14; 
options = odeset('RelTol',tol,'AbsTol',tol,'Refine',1);
tic; [tac,yac] = ode45(@fdpendulum,tspan,init,options); etac = toc;

% odepsmh
tic; [tpsm,ypsm] = odepsmh(@fdpendulum,tspan,stepsize,init,20); etpsm = toc;

% adaptive stepsizes
dtac = tac(2:end) - tac(1:end-1);  % dt = h's for accurate ode45
dtad = tad(2:end) - tad(1:end-1);  % dt of points for default ode45
dtad = 4*dtad(1:4:end); % actual since ode45 uses Refine 4 by default (see Help on odeset)

colormap(lines(4)); % attempting four colors distinguishable by colorblind
plot(yrk4(:,1),yrk4(:,2),...
     yad(:,1),yad(:,2),...
     yac(:,1),yac(:,2),...
     ypsm(:,1),ypsm(:,2))
legend(['rk4 h = ',num2str(stepsize)],'ode45 default','ode45 accurate',['odepsm h = ',num2str(stepsize)])
ax = gca;
ax.XTick = (-3:8)*pi;
ax.XTickLabel = cellstr([num2str((-3:8)'),repmat('\pi',12,1)]);
xlabel('y : clockwise angle from bottom');
ylabel('y'' : angular velocity');
grid on

% Find true (best) ending value using Symbolic Toolbox
% output from run below to avoid running this section of code
best = 17.417045282416465580661324578379; % with stepsize .25
%      17.417045282416465582192310843963 with stepsize .125
% using vpa (Variable Precision Arithmetic default 32 digits) 
% using odepsmh to degree 25 but it takes about 8 minutes on my old laptop
% stepsize2 = .25;
% tic; [tpsmvh,ypsmvh] = odepsmh(@fdpendulum,tspan,vpa(stepsize2),vpa(init),25); etpsmvh = toc;
% best = ypsmvh(end,1)  % assumed to be the most accurate value
% note: automatic conversion vpa(.1) from DE keeps the quadruple precision

% Table of results:
fprintf('%23s %10s %20s %8s\n','num steps','mean step','y(200) value','error');
str = ': %6d %10.4f %20.14f %8.1e\n'; 
fprintf(['  runge kutta 4',str],length(trk4)-1, stepsize, yrk4(end,1),abs(yrk4(end,1)-best));
fprintf(['  ode45 default',str],length(dtad), mean(dtad), yad(end,1),abs(yad(end,1)-best));
fprintf([' ode45 accurate',str],length(dtac), mean(dtac), yac(end,1),abs(yac(end,1)-best));
fprintf([' odepsmh deg 20',str],length(tpsm)-1, stepsize, ypsm(end,1),abs(ypsm(end,1)-best));