function dydt = ringModulator(t,y)
% Ring Modulator problemin "Test Set for IVP Solvers" Mazzi & Iavernaro 2003
% y(1:7) are voltages, y(8:15) are currents
% Constants
C = 1.6e-8;
Cs = 2e-12;
Cp = 1e-8;
Lh = 4.45;
Ls1 = 0.002;
Ls2 = 5e-4;
Ls3 = 5e-4;
gamma = 40.67286402e-9;
R = 25000;
Rp = 50;
Rg1 = 36.3;
Rg2 = 17.3;
Rg3 = 17.3;
Ri = 50;
Rc = 600;
delta = 17.7493332;

% reused expressions
Uin1 = 0.5*sin(2000*pi*t);
Uin2 =  2*sin(20000*pi*t);
UD1 = y(3) - y(5) - y(7) - Uin2;
UD2 = -y(4) + y(6) - y(7) - Uin2;
UD3 = y(4) + y(5) + y(7) + Uin2;
UD4 = -y(3) - y(6) + y(7) + Uin2;
qUD1 = gamma*(exp(delta*UD1)); % could overflow if delta*UDx > 172
qUD2 = gamma*(exp(delta*UD2));
qUD3 = gamma*(exp(delta*UD3));
qUD4 = gamma*(exp(delta*UD4));

dydt = y;  % initialize array of same class as y, overwrite values
dydt(1) = (y(8) - 0.5*y(10) + 0.5*y(11) + y(14) - y(1)/R)/C;
dydt(2) = (y(9) - 0.5*y(12) + 0.5*y(13) + y(15) - y(2)/R)/C;
dydt(3) = (y(10) - qUD1 + qUD4)/Cs;
dydt(4) = (-y(11) + qUD2 - qUD3)/Cs;
dydt(5) = (y(12) + qUD1 - qUD3)/Cs;
dydt(6) = (-y(13) - qUD2 + qUD4)/Cs;
dydt(7) = (-y(7)/Rp + qUD1 + qUD2 - qUD3 - qUD4)/Cp;
dydt(8) = -y(1)/Lh;
dydt(9) = -y(2)/Lh;
dydt(10) = (0.5*y(1) - y(3) - Rg2*y(10))/Ls2;
dydt(11) = (-0.5*y(1) + y(4) - Rg3*y(11))/Ls3;
dydt(12) = (0.5*y(2) - y(5) - Rg2*y(12))/Ls2;
dydt(13) = (-0.5*y(2) + y(6) - Rg3*y(13))/Ls3;
dydt(14) = (-y(1) + Uin1 - (Ri + Rg1)*y(14))/Ls1;
dydt(15) = (-y(2) - (Rc + Rg1)*y(15))/Ls1;