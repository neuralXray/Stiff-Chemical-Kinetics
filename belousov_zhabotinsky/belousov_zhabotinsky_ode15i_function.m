function dydt = belousov_zhabotinsky_ode15i_function(t, y, yp)
%% Chemical Parameters definition
% rate constants
k1 = 4.72;
k2 = 3e9;
k3 = 1.5e4;
k4 = 4e7;
k5 = 1;

dydt = [yp(1) + k1*y(1)*y(2),
        yp(2) + k1*y(1)*y(2) + k2*y(3)*y(2) - k5*y(6),
        yp(3) + k2*y(3)*y(2) - k3*y(3)*y(5) + 2*k4*y(3)^2 - k1*y(1)*y(2),
        yp(4) - k2*y(3)*y(2),
        yp(5) + k3*y(3)*y(5),
        yp(6) - k3*y(3)*y(5) + k5*y(6),
        yp(7) - k4*y(3)^2];

end

