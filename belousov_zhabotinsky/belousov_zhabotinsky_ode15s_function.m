function dydt = belousov_zhabotinsky_ode15s_function(t, y)
%% Chemical Parameters definition
% rate constants
k1 = 4.72;
k2 = 3e9;
k3 = 1.5e4;
k4 = 4e7;
k5 = 1;

dydt = [- k1*y(1)*y(2);
        - k1*y(1)*y(2) - k2*y(3)*y(2) + k5*y(6);
        - k2*y(3)*y(2) + k3*y(3)*y(5) - 2*k4*y(3)^2 + k1*y(1)*y(2);
        k2*y(3)*y(2);
        - k3*y(5)*y(3);
        k3*y(5)*y(3) - k5*y(6);
        k4*y(3)^2];

end

