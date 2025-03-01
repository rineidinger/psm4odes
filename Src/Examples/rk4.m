function [t,y] = rk4(fhandle,tspan,h,y0)
% Runge-Kutta 4 method for DE defined for one-column vector y'=f(t,y) where 
%   f is specified by any function name passed in as @funcname to fhandle,
%   tspan = [t0,tend],
%   h is the stepsize,
%   y0 is the initial column of y values at t0.
% To match ode45 format, f must return a column, and program returns tables
% where each row of t and y contain values at one step of the method.
% Richard Neidinger, 4/19/11, changed transpose from ' to .' 6/14/24  
t0 = tspan(1);
tend = tspan(2);
t = (t0:h:tend).'; % so t is a column vector of times t(k), where t(1)=t0.
if t(end) < tend
    t = [t;tend];
end
numpts = length(t);

y = zeros(numpts, length(y0));  % row y(k,:) will be the y-vector at t(k)
y(1,:) = y0.';
for k = 1:(numpts-1)
    if k == numpts-1; h = t(k+1)-t(k); end
    yk = y(k,:).'; % most recent row of y transposed to column vector
    q1 = feval(fhandle,t(k),yk);  % using old values of t and y
    q2 = feval(fhandle,t(k)+h/2,yk+h*q1/2);
    q3 = feval(fhandle,t(k)+h/2,yk+h*q2/2);
    q4 = feval(fhandle,t(k+1),yk+h*q3);
    newy = yk + (h/6)*(q1+2*q2+2*q3+q4);
    y(k+1,:) = newy.';  % transpose back into next row
end