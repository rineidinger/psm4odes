% Example using flame DE to show usage of psm4odes software
% Models flame from lighting a match. A famously stiff DE.
% https://www.mathworks.com/company/newsletters/articles/stiff-differential-equations.html
% Richard Neidinger 3/6/25
y0 = 1e-4; % 1e-n is progressively harder for n = 1,2,3,4,5,6,7
endt = 2/y0; % abrupt cliff around time 1/y0

deg = 20; % for one poly plot and odepsmh
h = 2; % stepsize for odepsmh (2 gives accuracy even as time span lengthens)
tol = 1e-14; % tolerance for odepsmJZ (1e-14 for dramatic improvement)
tolode = 2.22045e-14; % lowest RelTol allowed

makepsmcode(@fflame,1);  % writes or overwrites fflameseries.m

% The following four lines of code require and use Symbolic Toolbox
% showing how you can generate symbolic series coefficients
if license('test','Symbolic_Toolbox')
    degsym = 5;  % degree of symbolic results
    symcoefs = fflameseries(0,sym('a'),degsym);  % y0 = sym('a');
    symcoefs = expand(symcoefs).'
    vpacoefs = fflameseries(0,vpa(y0),degsym).'  % default quadruple precision
end

% one Maclaurin polynomial of degree deg
coefs = fflameseries(0,y0,deg);
ts = linspace(0,endt,100);
ys = serieseval(coefs,0,ts);  % poynomial values

% ODE Solvers
tic; [t45d,y45d] = ode45(@fflame,[0,endt],y0); et45d = toc; % default settings
t45d = t45d(1:3:end); % take out default Refine 4 interpolations
y45d = y45d(1:3:end);

options = odeset('RelTol',tolode,'AbsTol',tolode,'Refine',1);
tic; [t45,y45] = ode45(@fflame,[0,endt],y0,options); et45 = toc;
tic; [t15s,y15s] = ode15s(@fflame,[0,endt],y0,options); et15s = toc; % default settings% psm solvers
tic; [t23s,y23s] = ode23s(@fflame,[0,endt],y0, options); et23s = toc; 

% psm solvers
tic; [tpsmh,ypsmh] = odepsmh(@fflame,[0,endt],h,y0,deg); etpsmh = toc;
tic; [tpsm,ypsm,degJZ] = odepsmJZ(@fflame,[0,endt],y0,tol); etpsm = toc;

plot(t45d,y45d,t45,y45,ts,ys,tpsmh,ypsmh,'--',tpsm,ypsm,'.-',t15s,y15s,':',t23s,y23s,'-.');
axis([0,endt,-.1,1.1]);
legend('ode45 default',...
       ['ode45 tol = ',num2str(tolode)],...
       ['one poly deg = ',num2str(deg)],...
       ['odepsmh h = ',num2str(h),' deg = ',num2str(deg)],...
       ['odepsmJZ tol = ',num2str(tol)],...
       ['ode15s tol = ',num2str(tolode)],...
       ['ode23s tol = ',num2str(tolode)],...
       'location','Northwest');

% report results
str = ' used %7d steps in %5.3f seconds ';
fprintf(['   ode45',str,'with default options\n'], length(t45d)-1, et45d);
%fprintf(['  ode23s',str,'with default options\n'], length(t23s)-1, et23s);
fprintf(['   ode45',str,'for tol %7.1e\n'], length(t45)-1, et45, tolode);
fprintf(['  ode15s',str,'for tol %7.1e\n'], length(t15s)-1, et15s, tolode);
fprintf(['  ode23s',str,'for tol %7.1e\n'], length(t23s)-1, et23s, tolode);
fprintf([' odepsmh',str,'using degree %1d\n'], length(tpsmh)-1, etpsmh, deg);
fprintf(['odepsmJZ',str,'using degree %1d for tol %7.1e\n'], length(tpsm)-1, etpsm, degJZ, tol);
disp('Both ode45 and odepsmJZ perform best at very small tol when y0 is very small.')
disp('   Zoom around upper corner to see stiffness effect.')
