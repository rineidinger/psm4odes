function dqdt = pleiadescalar(t,q)
% Pleiades test problem (Mazzi & Iavernaro 2003)
% version using only scalar operations for PSM4ODES, compare vector version
% https://www.mathworks.com/help/matlab/math/solve-celestial-mechanics-problem-with-high-order-ode-solvers.html
  x = q(1:7);  % position of 7 stars
  y = q(8:14);
  dqdt = q; % preallocate same size and class as q
  dqdt(1:14) = q(15:28);  % dqdt(1:7) = x', dqdt(8:14) = y'.
  for i = 1:7
      xsum = 0;
      ysum = 0;      
      for j = [1:i-1,i+1:7]
          xDist = x(j)-x(i);
          yDist = y(j)-y(i);
          r = (xDist^2 + yDist^2)^(3/2);
          xsum = xsum + j*xDist/r;
          ysum = ysum + j*yDist/r;
      end
      dqdt(i+14) = xsum;  % dqdt(15:21) = x''
      dqdt(i+21) = ysum;  % dqdt(22:28) = y''
  end
end