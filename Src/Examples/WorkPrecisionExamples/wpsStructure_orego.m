% Comparing odepsmJZ to MATLAB solvers for a test problem.
% Creates work-precision diagram and steps-precision diagram.
% Data generated for many ode runs and stored in sturcture "solved".
% Richard Neidinger 12/19/25

% WARNING: This may take hours, so times just one run for each solver and
% tol.  Change once to false to use timeit to average run time over 
% multiple runs, but plan to run overnight.

titlename = 'OREGO Oregonator chemical model';
dename = 'orego';
dehandle = str2func(dename);
dim = 3; 
tspan = [0,360];
init = [1;2;3];
refsoln = [0.1000814870318523e1, 0.1228178521549917e4, 0.1320554942846706e3];  % from Test Set for IVP Solvers

solvers = ["odepsmJZ","ode45","ode113","ode15s","ode78","ode89"];
% assumes solvers(1) is "odepsmJZ", any number of others

tolerances = 2.3*10.^-(3:14); % 2.3 keeps above minimum relative tolerance 2.22045e-14 of MATLAB ode solvers
once = true;  % false uses timeit to average over multiple calls, which could take prohibative computation time.
              % true uses tic toc to time one run of each solver & tol, which takes much less computation time.

% Script saves solved data structure in filename ['wpdata_',dename,'.mat'].
% Recover solved with load filename.
% Rest of this script should remain same for other test problems.

clear solved % to start new structure
numtols = length(tolerances);
numsolvers = length(solvers);
for i = numsolvers:-1:1 % backwards preallocates on first time through
    solved(i).name = solvers(i);
    solved(i).tols = tolerances;
    solved(i).deg = NaN(1,numtols,'double'); % only needed for odepsmJZ
    solved(i).steps = NaN(1,numtols,'double');
    solved(i).scd = NaN(1,numtols,'double');
    solved(i).time = NaN(1,numtols,'double');
    solved(i).yend = NaN(numtols,length(refsoln),'double');
end

f4time = @() makepsmcode(dehandle,dim);
makecodetime = timeit(f4time);  % timeit averages across many runs
fprintf('time to generate psmcode = %10.4f\n',makecodetime)
format compact  % runode calls will display progress of multiple solves

for n = 1:numtols
    fprintf('Computing for tolerance %.1e\n',tolerances(n))
    options = odeset('RelTol',tolerances(n),'AbsTol',tolerances(n),'Refine',1);
    for i = 1:numsolvers
        [solved(i).yend(n,:), solved(i).scd(n), solved(i).time(n), solved(i).steps(n), solved(i).deg(n)]...
            = runode(str2func(solvers(i)),dehandle,tspan,init,options,refsoln,once);
    end
end

% return work precision steps data from all runs in a structure
save(['wpdata_',dename,'.mat'])

plotwps(solved, 1:numsolvers, titlename, makecodetime);  % may add arguments for manual axislimits
