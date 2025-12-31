% Comparing odepsmJZ to MATLAB solvers for a test problem.
% Creates work-precision diagram and steps-precision diagram.
% Data generated for many ode runs and stored in sturcture "solved".
% Uses timeit to average run time over multiple runs of same solver and tol
% Richard Neidinger 12/8/25

titlename = 'Pleiades Ceslestial Mechanics Problem';
dename = 'pleiadescalar';
dehandle = str2func(dename);
dim = 28; 
tspan = [0,3];  %"Test Set for IVP Solvers" is [0,2] in plot, [0,3] in ref
init = [3 3 -1 -3 2 -2 2 ...
        3 -3 2 0 0 -4 4 ...
        0 0 0 0 0 1.75 -1.5 ...
        0 0 0 -1.25 1 0 0]';
% scd (significant correct digits) based on the first 14 (position) components
refsoln = [  % from "Test Set for IVP Sovers"
  0.3706139143970502,...
  3.237284092057233,...
 -3.222559032418324,...
  0.6597091455775310,...
  0.3425581707156584,...
  1.562172101400631,...
 -0.7003092922212495,...
 -3.943437585517392,...
 -3.271380973972550,...
  5.225081843456543,...
 -2.590612434977470,...
  1.198213693392275,...
 -0.2429682344935824,...
  1.091449240428980
];
solvers = ["odepsmJZ","ode45","ode113","ode15s","ode78","ode89"];
% assumes solvers(1) is "odepsmJZ", any number of others

tolerances = 2.3*10.^-(3:14); % 2.3 keeps above minimum relative tolerance 2.22045e-14 of MATLAB ode solvers
once = false;  % false uses timeit to average over multiple calls, which could take prohibative computation time.
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

plotwps(solved, 1:numsolvers, titlename, makecodetime,[0,12.5,-.1,2.5],[0,12.5,-100,2500]);  % may add arguments for manual axislimits
