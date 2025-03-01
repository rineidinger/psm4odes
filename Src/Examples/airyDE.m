function dydt = airyDE(t,y)
% Airy's Equation:  y'' = t*y
% y(1) = y; y(2) = y':  y(1)' = y(2); y(2)' = t*y(1);
dydt = [ y(2); t*y(1) ];
