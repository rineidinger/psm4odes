% Example of generating symbolic series coefficients on fex1 DE
% also using vpa (Variable Precision Arithmetic)
% R. Neidinger, 12/6/24...12/24/25 (improved labeling)
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

y0 = 1/vpa(10);  % default quadruple precision is 32 digits
vpacoefs = fex1series(0,y0,5).'  % first 5 coefs to quadruple precision

echo off
% vpa works in odepsmh and odepsmJZ but it takes longer (about minute here)
disp('Computing odepsmh vpa solution on [0,12]...')
tspan = [0,12];
h = 1/vpa(8);  % h = .125 gives same result (converts to vpa) but take care
deg = 20;
[tpsmh,ypsmh] = odepsmh(@fex1,tspan,h,y0,deg);
fprintf('y(end) value is %s\n', ypsmh(end)) % use %s to print vpa digits
fprintf('using %g equal steps h=%f of degree %g\n\n', ...
    length(tpsmh)-1, h, deg)

disp('Computing odepsmJZ vpa solution on [0,12]...')
tol = 1e-20;
[tpsm,ypsm,degfortol] = odepsmJZ(@fex1,tspan,y0,tol);
fprintf('y(end) value is %s\n', ypsm(end)) %.32f would be double precision
fprintf('using tol %.1e and %g variable steps of degree %g\n', ...
    tol, length(tpsm)-1, degfortol)

% Both end values have about 19 significant correct digits (scd).  
% However, odepsmJZ does not improve for lower tolerances.  
% But odepsmh does improve by halving stepsize, quickly converging to 
% 1.7723906256380987828403493181175