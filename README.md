# PSM4ODES
by Richard D. Neidinger, Davidson College, 3/14/25.   
In MATLAB, generates series recurrence code for use in the Power Series Method (PSM) (or arbitrary-order Taylor series) to solve any Ordinary Differential Equations (ODEs) system of form y' = f(t,y) and y(t0) = y0.  
Accompanies preprint article *Automatic Series Recurrence Relations for Ordinary Differential Equations* and has a *userManual.pdf*.

# Contents and Usage:

PSM4ODES\Src\Functions folder contains the software codes and should be added to the MATLAB path by addpath command or Set Path environment button.   
Define a DE file, say f.m, as is usually done for any MATLAB ODE solver, where f takes and returns column vectors. For example, y'' = sin(y2) uses y1 = y and y2 = y' in the system y1' = y2 and y2' = sin(y1^2). Avoid preallocating output with zeros or explicit doubles. This f.m file  
&nbsp; &nbsp; &nbsp;   function dydt = f(t,y)  
&nbsp; &nbsp; &nbsp;   dydt = [y(2); sin(y(1)^2)];  
works well. Codes that start with dydt = zeros(2,1) will fail but may be made compatible by changing the preallocation to dydt = y, so that dydt replicates the class of y in the overloaded execution of f that generates the series recurrence code. The name f and names of the variables in f.m can be anything though t and y will be used in the generated series code. We expect and implement only scalar functions (so ^ not .^ in MATLAB).

ODE software
-----------------------
- `odepsmh.m`    Solves ode by series to order deg with step-size h,  
  * CALL:  `[t,Y] = odepsmh(@f,[t0,tend],h,Y0,deg);`

- `odepsmJZ.m`    Solves ode by series, tolerance determines deg and variable h,  
  * CALL:  `[t,Y] = odepsmJZ(@f,[t0,tend],Y0,tolerance);`

SERIES software (automatically called as needed by above)
---------------
- `makepsmcode.m`  Writes file `fseries` that generates coefs for ODE solution,  
  * CALL:  `makepsmcode(@f,numDEs)`  
  * CALL:  `coefs = fseries(t0,Y0,deg)`  where `f` is litearlly the name of the input `@f`.  
- `serieseval.m`  Evaluates series at ts using coefs about t0,  
  * CALL:  `values = serieseval(coefs,t0,ts)`

- `trace.m`  Class of objects used by makepsmcode:  
  If t and y are trace objects, each operation (arithmetic or transcendental) in evaluation of f(t,y) is overloaded (as defined in trace.m) to write the corresponding code for the evaluation and for the series recurrence relation.

Examples Folder Contents
--------------------------
- `airyDE.m`  defines system for Airy's equation.
- `airyDE_Script.m`  shows use of series coefficients and odepsmJZ.
  
- `fex1.m`  defines the first order differential equation.
- `fex1_Script.m`  shows use of series coefficients and odepsmJZ.
- `fex1_SymbolicScript.m`  computes symbolic and vpa coefficients and vpa solution.
  
- `fex2.m`  defines the system for the second order differential equation.
- `fex2_Script.m`  shows use of y and y′ series coefficients and odepsmJZ.

- `fdpendulum.m`  defines system for a Forced Damped Pendulum.
- `fdpendulum_Script.m`  to graphically compare odepsmh with ode45.
- `fdpendulum_Timing.m`  compares runtime for different solvers.

- `fflame.m`  defines the Flame Equation, a stiff first order differential equation.
- `fflame_Script.m`  compares solvers using low tolerance.
  
- `rk4.m`  classicRunge Kutta order 4 method for comparison.
