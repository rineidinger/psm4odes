% Example using flame DE to show usage of psm4odes software
% Models flame from lighting a match.
% https://www.mathworks.com/company/newsletters/articles/stiff-differential-equations.html
y0 = 1/30;
endt = 2/y0;
[t45,y45] = ode45(@fflame,[0,endt],y0);

makepsmcode(@fflame,1);  % writes or overwrites fflameseries.m
deg = 20;
coefs = fflameseries(0,y0,deg);
ts = 0:.5:endt;
ys = serieseval(coefs,0,ts);  % values of one polynomial of degree deg

h = 2; % stepsize
[tpsm,ypsm] = odepsmh(@fflame,[0,endt],h,y0,deg);

plot(t45,y45,ts,ys,tpsm,ypsm,'--');
axis([0,endt,0,1.1]);
legend('ode45 default',...
       ['one poly deg = ',num2str(deg)],...
       ['odepsm h = ',num2str(h)],  'location','Northwest');

% The following requires and uses Symbolic Toolbox
% showing how you can generate symbolic series coefficients
symcoefs = fflameseries(0,sym('a'),4);  % y0 = sym('a');
symcoefs = expand(symcoefs).'

vpacoefs = fflameseries(0,vpa(y0),4).'  % default quadruple precision