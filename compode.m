function [t45,y45,tpsm,ypsm,deg] = compode(fhandle,tspan,Y0,tol,plotvar)
% compare ode45 solver with odepsmJZ solver
% Richard Neidinger 2/2/21; 4/27/21 added Refine 1; 
% 10/14/22 error changed to max on interval vs just last point
% 10/25/22 allow psm that takes only one step, 7/19/24 simplified

options = odeset('RelTol',tol,'AbsTol',tol,'Refine',1);
tic; [t45,y45] = ode45(fhandle,tspan,Y0,options); et45 = toc;
tic; [tpsm,ypsm,deg] = odepsmJZ(fhandle,tspan,Y0,tol); etpsm = toc;

% only plots and compares the plotvar coordinate
if nargin < 5; plotvar = 1; end
plot(t45,y45(:,plotvar),tpsm,ypsm(:,plotvar),'--');
legend('ode45','odepsmJZ');

% Find good values to compare with for error
tol = 2.22045e-14;  % lowest relative tolerance allowed in ode45
options = odeset('RelTol',tol,'AbsTol',tol);
[t45,y45a] = ode45(fhandle,t45,Y0,options); % not reliable to tolerance
[tpsma,ypsma] = ode45(fhandle,tpsm,Y0,options);
if length(tpsm) == 2 % trival case of only one psm step
    ypsma = ypsma([1,end],:);
end

% report results
dt45 = t45(2:end) - t45(1:end-1); 
dtpsm = tpsm(2:end) - tpsm(1:end-1); 
fprintf('%16s %10s %9s %9s %8s\n','num steps','mean step','error','rel err','sec');
str = ': %6d %10.4f %9.2e %9.2e %8.3f\n';
y45err = max(abs(y45(:,plotvar)-y45a(:,plotvar)));
y45relerr = max(abs(y45(:,plotvar)-y45a(:,plotvar))./abs(y45a(:,plotvar)));
ypsmerr = max(abs(ypsm(:,plotvar)-ypsma(:,plotvar)));
ypsmrelerr = max( abs( (ypsm(:,plotvar)-ypsma(:,plotvar))./ypsma(:,plotvar) ) );
fprintf(['   ode45',str],length(dt45), mean(dt45), y45err, y45relerr, et45);
fprintf(['odepsmJZ',str],length(dtpsm), mean(dtpsm), ypsmerr, ypsmrelerr, etpsm);
fprintf('odepsmJZ used series of degree %d\n',deg);