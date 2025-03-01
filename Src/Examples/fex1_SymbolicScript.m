% Example of generating symbolic series coefficients on fex1 DE
% also using vpa (Variable Precision Arithmetic)
% R. Neidinger, 12/6/24...2/21/25.
% Assumes fex1series.m has already been created by makepsmcode(@fex1,1).
% The following requires and uses MATLAB's Symbolic Toolbox.
echo on
y0 = sym('a');
symcoefs = fex1series(0,y0,3)
symcoefs.'  % transpose may help in reading output
% expand and simplify commands may also help

y0 = 1/sym(10);  % exact numerical
% Note: sym(1/10) works but 100000/sym(1234567) requires this better form
symcoefs = fex1series(0,y0,5)  % first 5 coefs exact numerical
disp(string(symcoefs.'))

y0 = 1/vpa(10);  % default quadruple precision
vpacoefs = fex1series(0,y0,5).'  % first 5 coefs to quadruple precision

echo off
% vpa works in odepsmh and odepsmJZ but it takes longer (about minute here)
disp('Computing odepsmh solution using variable precision...')
tspan = [0,12];
deg = 20;
y0 = 1/vpa(10); % default precision is 32 digits
h = 1/vpa(8);  % h = .125 gives same result (converts to vpa) but take care
[tpsmh,ypsmh] = odepsmh(@fex1,tspan,h,y0,deg);
disp('y(end) value is')
better = ypsmh(end)
fprintf(['odepsmh used ', ...
    num2str(length(tpsmh)-1),' equal steps of deg ', num2str(deg),'.\n'])

disp('Computing odepsmJZ solution using variable precision...')
tol = 1e-24;
[tpsm,ypsm,degfortol] = odepsmJZ(@fex1,tspan,y0,tol);
best = ypsm(end)
fprintf(['For tol ',num2str(tol),', odepsmJZ used ', ...
    num2str(length(tpsm)-1),' variable steps of deg ', num2str(degfortol),'.\n'])
