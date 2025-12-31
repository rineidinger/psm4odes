function dydt = orego(t,y)
% OREGO problem from Test Set for IVP Problems, Mazzia & Iavernaro 2003/08
s = 77.27;
w = 0.161;
q = 8.375e-6;
dydt = y; % preallocate y with same class as y
dydt(1) = s*(y(2) - y(1)*y(2) + y(1) - q*y(1)^2);
dydt(2) = (-y(2) - y(1)*y(2) + y(3))/s;
dydt(3) = w*(y(1) - y(3));