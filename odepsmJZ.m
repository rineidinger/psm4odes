function [t,y,deg] = odepsmJZ(fhandle,tspan,y0,tolerance)
% Adaptive power series method for solving scalar DE y'=f(t,y) with
% variable stepsize and fixed order determined by tolerance.
%   f is specified by any function name passed in as @funcname to fhandle,
%   tspan = [t0,tend],
%   y0 is the initial y values at t0, vector column to match f in and out,
%   tolerance is the desired absolute and relative error bound.
%
% Uses order and step-size estimation from [JZ] A. Jorba & M. Zou, 
% A Software % Package for the Nuerical Integration of ODEs by Means of 
% High-Order Taylor Methods, Experimental Mathematics, 14:1, 99-117.
%
% Uses makepsmcode.m and trace.m which must be in path.
% Writes and calls psm file [func2str(fhandle),'series.m'] to generate series.
%
% To match ode45 format, f must return a column, and program returns tables
% where each row of t and y contain values at one step of the method.
% Internally, y components will go down columns as used in f, before final
% transpose for output.
% Richard Neidinger 2/2/21; update 5/4/21 for min step and time direction
t0 = tspan(1);
tend = tspan(2);
numDEs = length(y0);  % get number of DEs in f from y0
y0 = reshape(y0,[numDEs,1]);  % insure a column vector
t = t0; % row vector that will grow by variable stepsize
y = y0;  % y0 is first column, will augment more columns for each t
hsign = sign(tend - t0);  % direction of time
minstep = 2^-10; % minimum stepsize to prevent infinite loop
minstepflag = false;

seriesfuncname = [func2str(fhandle),'series'];
if isfile([seriesfuncname,'.m'])
    fseries = str2func(seriesfuncname);  % return handle to that function
else
    fseries = makepsmcode(fhandle,numDEs);
end

% fixed order
deg = ceil(1 - log(tolerance)/2); % [JZ] (3-4) w/ eps_a = eps_r = tolerance

k = 1;
while hsign* (t(k) - tend) < 0
    coefs = fseries(t(k),y(:,k),deg);
    
    % determine stepsize by [JZ] root test algorithm
    scale = max(1,norm(y(:,k),inf));
    roottests = (scale./vecnorm(coefs(:,2:end),inf,1)).^(1./(1:deg));
    hm = min(roottests([deg-1,deg]))/exp(2) * exp(-0.7/(deg-1));
    h = min(hm, min(roottests));
    if h < minstep
        h = minstep;
        if ~minstepflag
            note = [num2str(minstep),', first at t = ', num2str(t(k))];
            disp(['Warning: odepsmJZ solver requested below minimum stepsize ',note]);
            minstepflag = true;
        end
    end
    h = min(h, abs(tend - t(k))); % only so last step hits tend
    t(k+1) = t(k) + hsign*h;  % augment t vector
    
    % horners rule to evaluate polynomial in powers of (t-t(k)) at t(k+1)
    val = coefs(:,deg+1);
    for m = deg:-1:1
        val = val*hsign*h + coefs(:,m);
    end
    y(:,k+1) = val;  % augment y array
    k = k+1;
end
y = y.'; % transpose into output form
t = t.';