function dydt = fflame(t,y)
% forward flame DE
dydt = y^2*(1 - y);