function dydt = hires(t,y)
% HIRES problem from Test Set for IVP Problems, Mazzia & Iavernaro 2003/08
dydt = y; % preallocate y with same class as y
dydt(1) = -1.71*y(1) + 0.43*y(2) + 8.32*y(3) + 0.0007;
dydt(2) = 1.71*y(1) - 8.75*y(2);
dydt(3) = -10.03*y(3) + 0.43*y(4) + 0.035*y(5);
dydt(4) = 8.32*y(2) + 1.71*y(3) - 1.12*y(4);
dydt(5) = -1.745*y(5) + 0.43*y(6) + 0.43*y(7);
dydt(6) = -280*y(6)*y(8) + 0.69*y(4) + 1.71*y(5) - 0.43*y(6) + 0.69*y(7);
dydt(7) = 280*y(6)*y(8) - 1.81*y(7);
dydt(8) = -280*y(6)*y(8) + 1.81*y(7);