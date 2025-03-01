% Series Method on fex1.m: dydt = sin(y^2)
% Richard Neidinger 12/5/24
t0 = 0; tend = 12;
tspan = [t0,tend];
y0 = 0.1;
tol = 1e-9;
deg = 25; % below assumes > 10

makepsmcode(@fex1,1);  % writes fex1series.m code 
% fex1series generates coefficients about point (t0,y0) to degreee deg
coefs = fex1series(t0,y0,deg);
ts = linspace(t0,tend,100);
ys = serieseval(coefs,t0,ts);  % values of degree deg Taylor poly solution
y5 = serieseval(coefs(1:6),t0,ts);  % values of degree 5 Taylor poly 
y10 = serieseval(coefs(1:11),t0,ts);  % values of degree 10 Taylor poly

options = odeset('RelTol',tol,'AbsTol',tol,'Refine',1); % 1 pt per step
[t45,y45] = ode45(@fex1,tspan,y0,options);

% adaptive piecewise series solution using fex1series:
[tpsm,ypsm,degfortol] = odepsmJZ(@fex1,tspan,y0,tol); 

plot(t45,y45,"k",ts,y5,"k--",ts,y10,"k:",ts,ys,"k-.",tpsm,ypsm,'k*');
legend('ode45','deg 5 Maclaurin poly','deg 10',['deg ',num2str(deg)],...
    ['piecewise deg ',num2str(degfortol)],'Location','northwest')
xlabel('t');
ylabel('y');
axis([t0,tend,0,2])  % y range hard coded, change for other parameters
                     % Maclaurin label hard coded for t0 = 0
fprintf(['For tol ',num2str(tol),', ode45 used ',num2str(length(t45)-1),...
    ' steps, odepsmJZ used ',num2str(length(tpsm)-1),' steps of deg ', ...
    num2str(degfortol),'.\n'])