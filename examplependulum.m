% Example comparing odepsmh and odepsmJZ to ode45 for forced damped pendulum
% Richard Neidinger, 1/18/2021, 4/27/21 Refine 1 option, 6/19/24 formatting
format compact
format shorte
echo on

stepsize = .6;  % et stands for elapsed time in sec
tic; [trk4,yrk4] = rk4(@fdpendulum,[0,200],stepsize,[0;2]); etrk = toc; 
tic; [tad,yad] = ode45(@fdpendulum,[0,200],[0;2]); etad = toc;
% default ode45 uses 'Refine' 4 that interpolates four values for each step
tic; [tpsm,ypsm] = odepsmh(@fdpendulum,[0,200],stepsize,[0;2],20); etpsm = toc;
options = odeset('RelTol',2.3e-14,'AbsTol',1e-15,'Refine',1);
tic; [tac,yac] = ode45(@fdpendulum,[0,200],[0;2],options); etac = toc;
dtac = tac(2:end) - tac(1:end-1);  % dt = h's for accurate ode45
dtad = tad(2:end) - tad(1:end-1);  % dt = h's for default ode45
dtad4 = 4*dtad(1:4:end); % actual since ode45 uses Refine 4 by default (see Help on odeset)

% try odepsmJZ with tolerance:
tol = 1e-13;
tic; [tpsmad,ypsmad,degad] = odepsmJZ(@fdpendulum,[0,200],[0;2],tol); etpsmad = toc;
dtpsmad = tpsmad(2:end) - tpsmad(1:end-1);  

stepsize2 = .25;
tic; [tpsm2,ypsm2] = odepsmh(@fdpendulum,[0,200],stepsize2,[0;2],20); etpsm2 = toc;
best = ypsm2(end,1);  % assumed to be the most accurate value
echo off

plot(yrk4(:,1),yrk4(:,2),'g',...
     yad(:,1),yad(:,2),'b',...
     yac(:,1),yac(:,2),'k',...
     ypsm(:,1),ypsm(:,2),'r')
legend(['rk4 h = ',num2str(stepsize)],'ode45 default','ode45 accurate',['odepsm h = ',num2str(stepsize)])
ax = gca;
ax.XTick = (-3:8)*pi;
ax.XTickLabel = cellstr([num2str((-3:8)'),repmat('\pi',12,1)]);
xlabel('y : clockwise angle from bottom');
ylabel('y'' : angular velocity');
grid on

fprintf('%23s %10s %20s %8s %8s\n','num steps','mean step','y(200) value','error','sec');
str = ': %6d %10.4f %20.14f %8.1e %8.1e\n';
fprintf(['  runge kutta 4',str],length(trk4)-1, stepsize, yrk4(end,1),abs(yrk4(end,1)-best),etrk);
fprintf(['  ode45 default',str],length(dtad4), mean(dtad4), yad(end,1),abs(yad(end,1)-best),etad);
fprintf([' ode45 accurate',str],length(dtac), mean(dtac), yac(end,1),abs(yac(end,1)-best),etac);
fprintf([' odepsmh deg 20',str],length(tpsm)-1, stepsize, ypsm(end,1),abs(ypsm(end,1)-best),etpsm);
fprintf([' odepsmh deg 20',str],length(tpsm2)-1, stepsize2, ypsm2(end,1),0,etpsm2);
fprintf(['odepsmJZ deg %2d',str],degad,length(dtpsmad), mean(dtpsmad), ypsmad(end,1),abs(ypsmad(end,1)-best),etpsmad);
% Although zero error is reported for this "best" value, actual error is
% unknown and relative error is not guaranteed to be within machine eps.
% try odepsmJZ with tolerance:
tol = 1e-3;
tic; [tpsmad,ypsmad,degad] = odepsmJZ(@fdpendulum,[0,200],[0;2],tol); etpsmad = toc;
dtpsmad = tpsmad(2:end) - tpsmad(1:end-1);
fprintf(['odepsmJZ deg %2d',str],degad,length(dtpsmad), mean(dtpsmad), ypsmad(end,1),abs(ypsmad(end,1)-best),etpsmad);