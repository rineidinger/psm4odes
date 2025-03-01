function dydt = fdpendulum(t,y)
% forced damped pendulum p'' = -sin(p) - .1*p' + cos(t)
% y(1) = p;  y(2) = p'
dydt = [ y(2); -sin(y(1)) - 0.1*y(2) + cos(t)];