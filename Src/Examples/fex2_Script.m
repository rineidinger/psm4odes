% Series Method on fex2.m: dydt = [y(2); sin(y(1)^2)];
% Richard Neidinger 12/11/24... 2/26/25
t0 = 0; tend = 12;
tspan = [t0,tend];
y0 = [1;0];  % [1,0] will also work
tol = 1e-9;
deg = 25; % below assumes > 10

makepsmcode(@fex2,2);  % writes fex2series.m code 
% fex2series generates coefficients about point (t0,y0) to degreee deg
coefs = fex2series(t0,y0,deg);  % two rows of coefs for y and y'

ts = linspace(t0,tend,100);
% 2 rows of values for y and y' at each ts
ys = serieseval(coefs,t0,ts);  % values of degree deg Taylor poly solution
y5 = serieseval(coefs(:,1:6),t0,ts);  % values of degree 5 Taylor poly 
y10 = serieseval(coefs(:,1:11),t0,ts);  % values of degree 10 Taylor poly
ys = ys.'; y5 = y5.'; y10 = y10.'; % transpose into 2 columns as ode output

options = odeset('RelTol',tol,'AbsTol',tol,'Refine',1); % 1 pt per step
[t45,y45] = ode45(@fex2,tspan,y0,options);

% adaptive piecewise series solution using fex2series:
[tpsm,ypsm,degfortol] = odepsmJZ(@fex2,tspan,y0,tol); 

% phase plane (y,y') plot
plot(y45(:,1),y45(:,2),"k",ypsm(:,1),ypsm(:,2),'k*');
hold on
axis manual  % keep scale of adaptive plots
% plot one Maclaurin series
plot(y5(:,1),y5(:,2),"g--",y10(:,1),y10(:,2),"b:",ys(:,1),ys(:,2),"r-.");
legend('ode45',['piecewise deg ',num2str(degfortol)], ...
    'deg 5 Maclaurin poly','deg 10',['deg ',num2str(deg)], ...
    'Location','southeast')
xlabel('y');
ylabel('y''');

fprintf(['For tol ',num2str(tol),', ode45 used ',num2str(length(t45)-1),...
    ' steps, odepsmJZ used ',num2str(length(tpsm)-1),' steps of deg ', ...
    num2str(degfortol),'.\n'])