% Comparing odepsmJZ to MATLAB solvers for a test problem.
% Creates work-precision diagram and steps-precision diagram.
% Data generated for many ode runs and stored in sturcture "solved".
% Richard Neidinger 12/19/25

% WARNING: This intense problem may take many hours to run.  We skip the
% lowest tolerance 2.3e-14 for ode78 and ode89, since they fail giving
% error exceeds maximum array size 31.7GB (using almost 284 million steps).
% This times only one run of each solver and tol since this is very time
% consuming with tens of millions of steps.

titlename = 'Ring Modulator';
dename = 'ringModulator';
dehandle = str2func(dename);
dim = 15; 
tspan = [0,1e-3];
init = zeros(15,1);
refsolntestset = [  % from Test Set for IVP Solvers
    -0.2339057358486745e-1,...
    -0.7367485485540825e-2,...
     0.2582956709291169,...
    -0.4064465721283450,...
    -0.4039455665149794,...
     0.2607966765422943,...
     0.1106761861269975,...
     0.2939904342435596e-6,...
    -0.2840029933642329e-7,...
     0.7267198267264553e-3,...
     0.7929487196960840e-3,...
    -0.7255283495698965e-3,...
    -0.7941401968526521e-3,...
     0.7088495416976114e-4,...
     0.2390059075236570e-4
];

solvers = ["odepsmJZ","ode45","ode113","ode15s","ode78","ode89"];
% assumes solvers(1) is "odepsmJZ", any number of others

tolerances = 2.3*10.^-(4:2:14); % 2.3 keeps above minimum relative tolerance 2.22045e-14 of MATLAB ode solvers
% Uses fewer tolerances than other test problems to save run time.
once = true;  % false uses timeit to average over multiple calls, which could take prohibative computation time.
              % true uses tic toc to time one run of each solver & tol, which takes much less computation time.

% Script saves solved data structure in filename ['wpdata_',dename,'.mat'].
% Recover solved with load filename.
% Rest of this script should remain same for other test problems

% For Ring Modulator, added line 68 if, to avoid error that stops program

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
        if ~(n == numtols && i > 4) % avoid ode78 and ode89 that exceed maximum array for 2.3e-14 tol
        [solved(i).yend(n,:), solved(i).scd(n), solved(i).time(n), solved(i).steps(n), solved(i).deg(n)]...
            = runode(str2func(solvers(i)),dehandle,tspan,init,options,refsoln,once);
        end
    end
end

% return work precision steps data from all runs in a structure
save(['wpdata_',dename,'.mat'])

plotwps(solved, 1:numsolvers, titlename, makecodetime);  % may add arguments for manual axislimits
