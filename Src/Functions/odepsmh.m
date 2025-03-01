function [t,y] = odepsmh(fhandle,tspan,h,y0,deg)
% power series method for solving scalar DE y'=f(t,y) w/ fixed stepsize 
%   f is specified by any function name passed in as @funcname to fhandle,
%   tspan = [t0,tend],
%   h is the fixed stepsize,
%   y0 is the initial y values at t0, vector column to match f in and out
%   deg is the degree of the Taylor polynomials to be computed.
%
% Uses makepsmcode.m and trace.m which must be in path.
% Writes and calls psm file [func2str(fhandle),'series.m'] to generate series.
%
% To match ode45 format, f must return a column, and program returns tables
% where each row of t and y contain values at one step of the method.
% Internally, y components will go down columns as used in f, before final
% transpose for output.
% Richard Neidinger 12/3/20...12/6/24 
t0 = tspan(1);
tend = tspan(2);
if h*(tend - t0) < 0; h = -h; end  % match sign of h to direction
t = (t0:h:tend).'; % t is a column vector of times t(k), with t(1)=t0
if t(end) ~= tend
    t = [t;tend];
end
numpts = length(t);

% series recurrence code is generated if it doesn't already exist
seriesfuncname = [func2str(fhandle),'series'];
if isfile([seriesfuncname,'.m'])
    fseries = str2func(seriesfuncname);  % return function handle
else
    fseries = makepsmcode(fhandle,length(y0));
end

% preallocate y with copies of y0; only first column values are used
% using y0 instead of zeros allows vpa y0 using Symbolic Toolbox
y = repmat(y0,1,numpts);
% indexed by y(component, time)

for k = 1:(numpts-1)
    coefs = fseries(t(k),y(:,k),deg);  
    % horners rule to evaluate polynomial in powers of (t(k+1)-t(k))
    val = coefs(:,deg+1);
    h = t(k+1)-t(k); % could differ on last step
    for m = deg:-1:1
        val = val*h + coefs(:,m);
    end
    y(:,k+1) = val;
end
y = y.'; % transpose into output form