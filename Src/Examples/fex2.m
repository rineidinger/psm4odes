function dydt = fex2(t,y)
% system for y'' = sin(y^2), so y(1) = y, y(2) = y'
dydt = [y(2); sin(y(1)^2)];