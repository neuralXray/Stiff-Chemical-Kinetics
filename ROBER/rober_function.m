function dydt = rober_function(y)
%% Chemical Parameters definition
% rate constants
k1 = 4e-2;
k2 = 3e7;
k3 = 1e4;

dydt = [- k1*y(:, 1) + k3*y(:, 2).*y(:, 3), ...
        k1*y(:, 1) - k2*y(:, 2).^2 - k3*y(:, 2).*y(:, 3), ...
        k2*y(:, 2).^2];

end

