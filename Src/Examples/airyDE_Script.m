% Series Method on airyDE.m: dydt = [ y(2); t*y(1) ];
% Richard Neidinger 2/27/25
t0 = 0; tend = -10; % going in negative direction is more interesting plot
tspan = [t0,tend];
y0 = [1;0];  % [1;0] and [0;1] give the two linearly independent solutions
tol = 1e-9;
deg = 25; % below assumes > 10

makepsmcode(@airyDE,2);  % writes airyDEseries.m code 
% airyDEseries generates coefficients about point (t0,y0) to degreee deg
coefs = airyDEseries(t0,y0,deg);  % two rows of coefs for y and y'
format compact
format shortE
disp('y coefs for specified initial conditions:')
disp(coefs(1,:))

ts = linspace(t0,tend,100);
% 2 rows of values for y and y' at each ts
ys = serieseval(coefs,t0,ts);  % values of degree deg Taylor poly solution
y5 = serieseval(coefs(:,1:6),t0,ts);  % values of degree 5 Taylor poly 
y10 = serieseval(coefs(:,1:11),t0,ts);  % values of degree 10 Taylor poly
ys = ys.'; y5 = y5.'; y10 = y10.'; % transpose into 2 columns as ode output

options = odeset('RelTol',tol,'AbsTol',tol,'Refine',1); % 1 pt per step
[t45,y45] = ode45(@airyDE,tspan,y0,options);

% adaptive piecewise series solution using airyDEseries:
[tpsm,ypsm,degfortol] = odepsmJZ(@airyDE,tspan,y0,tol); 

% plot y vs t using first column of y matrices
plot(t45,y45(:,1),"k",tpsm,ypsm(:,1),'k*');
hold on
axis manual  % keep scale of adaptive plots
% plot one Maclaurin series
plot(ts,y5(:,1),"g--",ts,y10(:,1),"b:",ts,ys(:,1),"r-.");
legend('ode45',['piecewise deg ',num2str(degfortol)], ...
    'deg 5 Maclaurin poly','deg 10',['deg ',num2str(deg)], ...
    'Location','southwest')
xlabel('t');
ylabel('y');

fprintf(['For tol ',num2str(tol),', ode45 used ',num2str(length(t45)-1),...
    ' steps, odepsmJZ used ',num2str(length(tpsm)-1),' steps of deg ', ...
    num2str(degfortol),'.\n'])

% The following four lines of code require and use Symbolic Toolbox
% Comment out these four lines if Symbolic Toolbox is not installed.
degsym = 15;  % degree of symbolic results
symcoefs = airyDEseries(0,[sym('a'),sym('b')],degsym);
disp('y coefs for y(0)=[a;b]:')
disp(symcoefs(1,:))