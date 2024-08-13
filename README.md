# psm4odes
by Richard D. Neidinger, Davidson College, 7/19/24.   
In MATLAB, generates series recurrence code for use in the Power Series Method (PSM) (or arbitrary-order Taylor series) to solve any Ordinary Differential Equations (ODEs) system of form y' = f(t,y) and y(t0) = y0.  
Accompanies article *Automatic Series Recurrence Relations for Ordinary Differential Equations* to be submitted to SIAM Review in 2024.

# Contents and Usage:

ODE is specified, just as for MATLAB ode45, by a function such as `f.m`, two such provided examples are:
- `fflame.m`
- `fdpendulum.m`

WARNING: DO NOT preallocate the array (vector) output in the ODE function using zeros or explicit doubles, instead construct it by concatenating variables or expressions or just preallocate with input y (all of whick can be overloaded with a trace object type).

ODE software
-----------------------
- `odepsmh.m`    Solves ode by series to order deg with step-size h,  
  * CALL:  `[t,Y] = odepsmh(@f,[t0,tend],h,Y0,deg);`

- `odepsmJZ.m`    Solves ode by series, tolerance determines deg and variable h,  
  * CALL:  `[t,Y] = odepsmJZ(@f,[t0,tend],Y0,tolerance);`

- `compode.m`    Solves ode by ode45 and odepsmJZ and compares results,
  * CALL:  `[t45,y45,tpsm,ypsm] = compode(@f,[t0,tend],Y0,tolerance);`

SERIES software (automatically called as needed by above)
---------------
- `makepsmcode.m`  Writes file `fseries` that generates coefs for ODE solution,  
  * CALL:  `makepsmcode(@f,numDEs)`  
  * CALL:  `coefs = fseries(t0,Y0,deg)`  where `f` is litearlly the name of the input `@f`.  
- `serieseval.m`  Evaluates series at ts using coefs about t0,  
  * CALL:  `values = serieseval(coefs,t0,ts)`

- `trace.m`  Class of objects used by makepsmcode:  
  If t and y are trace objects, each operation (arithmetic or transcendental) in evaluation of f(t,y) is overloaded (as defined in trace.m) to write the corresponding code for the evaluation and for the series recurrence relation.

Documentation and Examples
--------------------------
- `README.md`  This text file.

- `AutoSerRecur.pdf`  Preprint article to accompany this repository.

- `examplependulum.m`  Script to compare odepsmh with ode45 on fdpendulum.
  * CALL:  `examplependulum`
- `exampleflame.m`  Script showing direct use of psm code (incl. symbolical)
  * CALL:  `exampleflame`
- `rk4.m`  Runge Kutta order 4 method for comparison.
