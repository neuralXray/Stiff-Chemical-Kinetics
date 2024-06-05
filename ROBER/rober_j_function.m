function J = rober_j_function(h_0, hd, ci, yi)
%% Chemical Parameters definition
% rate constants
k1 = 4e-2;
k2 = 3e7;
k3 = 1e4;

Loss1_betai1 = ci*hd + k1*h_0;
Loss1_betai2 = - k3*yi(:, 3).*h_0;
Loss1_betai3 = - k3*yi(:, 2).*h_0;

Loss2_betai1 = - k1*h_0;
Loss2_betai2 = ci*hd + (2*k2*yi(:, 2) + k3*yi(:, 3)).*h_0;
Loss2_betai3 = k3*yi(:, 2).*h_0;

Z = zeros(size(h_0));
Loss3_betai2 = - 2*k2*yi(:, 2).*h_0;
Loss3_betai3 = ci*hd;

J = [Loss1_betai1, Loss1_betai2, Loss1_betai3;
     Loss2_betai1, Loss2_betai2, Loss2_betai3;
         Z       , Loss3_betai2, Loss3_betai3];

end

