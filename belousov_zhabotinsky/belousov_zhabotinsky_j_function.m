function J = belousov_zhabotinsky_j_function(h_0, hd, ci, yi)
%% Chemical Parameters definition
% rate constants
k1 = 4.72;
k2 = 3e9;
k3 = 1.5e4;
k4 = 4e7;
k5 = 1;

Loss1_betai1 = ci*hd + k1*yi(:, 2).*h_0;
Loss1_betai2 = k1*yi(:, 1).*h_0;
Z = zeros(size(h_0));

Loss2_betai1 = k1*yi(:, 2).*h_0;
Loss2_betai2 = ci*hd + k1*yi(:, 1).*h_0 + k2*yi(:, 3).*h_0;
Loss2_betai3 = k2*yi(:, 2).*h_0;       
Loss2_betai6 = - k5*h_0;

Loss3_betai1 = - k1*yi(:, 2).*h_0;
Loss3_betai2 = k2*yi(:, 3).*h_0 - k1*yi(:, 1).*h_0; 
Loss3_betai3 = ci*hd + k2*yi(:, 2).*h_0 - k3*yi(:, 5).*h_0 + 4*k4*yi(:, 3).*h_0;        
Loss3_betai5 =  -k3*yi(:, 3).*h_0;

Loss4_betai2 = - k2*yi(:, 3).*h_0;
Loss4_betai3 = - k2*yi(:, 2).*h_0;
Loss4_betai4 = ci*hd;

Loss5_betai3 = k3*yi(:, 5).*h_0;
Loss5_betai5 = ci*hd + k3*yi(:, 3).*h_0;

Loss6_betai3 = - k3*yi(:, 5).*h_0;
Loss6_betai5 = - k3*yi(:, 3).*h_0;
Loss6_betai6 = ci*hd + k5*h_0;

Loss7_betai3 = - 2*k4*yi(:, 3).*h_0;
Loss7_betai7 = ci*hd;

% Jacobian matrix     
J = [Loss1_betai1, Loss1_betai2,      Z      ,      Z      ,      Z      ,      Z      ,      Z      ; 
     Loss2_betai1, Loss2_betai2, Loss2_betai3,      Z      ,      Z      , Loss2_betai6,      Z      ;
     Loss3_betai1, Loss3_betai2, Loss3_betai3,      Z      , Loss3_betai5,      Z      ,      Z      ;
          Z      , Loss4_betai2, Loss4_betai3, Loss4_betai4,      Z      ,      Z      ,      Z      ;
          Z      ,      Z      , Loss5_betai3,      Z      , Loss5_betai5,      Z      ,      Z      ;
          Z      ,      Z      , Loss6_betai3,      Z      , Loss6_betai5, Loss6_betai6,      Z      ;
          Z      ,      Z      , Loss7_betai3,      Z      ,      Z      ,      Z      , Loss7_betai7];

end

