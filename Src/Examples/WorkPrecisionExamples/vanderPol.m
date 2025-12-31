function dydt = vanderPol(t,y)
% van der Pol RLC oscillator z'' = mu*(1-z^2)*z' - z
% y(1) = z;  y(2) = z'
mu = 1000;
dydt = [ y(2); mu*(1-y(1)^2)*y(2) - y(1) ];