function J = akzonobel_j_function(h_0, hd, ci, yi)
%% Chemical Parameters definition
% rate constants
k1 = 18.7;
k2 = 0.58;
k3 = 0.09;
k4 = 0.42;
K = 34.4;
klA = 3.3;
pCO2 = 0.9;
H = 737;
        
% L1
Loss1_betai1 = ci*hd + (8*k1*yi(:, 1).^3.*yi(:, 2).^0.5 + (k2/K)*yi(:, 5) + k3*yi(:, 4).^2).*h_0;
Loss1_betai2 = 0.5*k1*yi(:, 1).^4.*yi(:, 2).^-0.5.*h_0;
Loss1_betai3 = -k2*yi(:, 4).*h_0;
Loss1_betai4 = (-k2*yi(:, 3) + 2*k3*yi(:, 1).*yi(:, 4)).*h_0;
Loss1_betai5 = (k2/K)*yi(:, 1).*h_0;
Z = zeros(size(h_0));

%L2
Loss2_betai1 = (2*k1*yi(:, 1).^3.*yi(:, 2).^0.5 + k3*yi(:, 4).^2).*h_0;
Loss2_betai2 = ci*hd + (0.25*k1*yi(:, 1).^4.*yi(:, 2).^-0.5 + 0.25*k4*yi(:, 6).^2.*yi(:, 2).^-0.5 + klA).*h_0;      
Loss2_betai4 = 2*k3*yi(:, 1).*yi(:, 4).*h_0;
Loss2_betai6 = k4*yi(:, 6).*yi(:, 2).^0.5.*h_0;

%L3
Loss3_betai1 = (-4*k1*yi(:, 1).^3.*yi(:, 2).^0.5 - (k2/K)*yi(:, 5)).*h_0;
Loss3_betai2 = -0.5*k1*yi(:, 1).^4.*yi(:, 2).^-0.5.*h_0; 
Loss3_betai3 = ci*hd + k2*yi(:, 4).*h_0;        
Loss3_betai4 = k2*yi(:, 3).*h_0;
Loss3_betai5 = - (k2/K)*yi(:, 1).*h_0;

%L4
Loss4_betai1 = ( -(k2/K)*yi(:, 5) + 2*k3*yi(:, 4).^2).*h_0;
Loss4_betai3 = k2*yi(:, 4).*h_0;
Loss4_betai4 = ci*hd + (k2*yi(:, 3) + 4*k3*yi(:, 1).*yi(:, 4)).*h_0;
Loss4_betai5 = (k2/K)*yi(:, 1).*h_0;

%L5
Loss5_betai1 = (k2/K)*yi(:, 5).*h_0;
Loss5_betai2 = -0.5*k4*yi(:, 6).^2.*yi(:, 2).^-0.5.*h_0;
Loss5_betai3 = -k2*yi(:, 4).*h_0;
Loss5_betai4 = -k2*yi(:, 3).*h_0;
Loss5_betai5 = ci*hd + (k2/K)*yi(:, 1).*h_0;
Loss5_betai6 = -2*k4*yi(:, 6).*yi(:, 2).^0.5.*h_0;

%L6
Loss6_betai2 = 0.5*k4*yi(:, 6).^2.*yi(:, 2).^-0.5.*h_0;
Loss6_betai6 = ci*hd + 2*yi(:, 6).*yi(:, 2).^0.5.*h_0;
       
% Jacobian matrix     
J = [Loss1_betai1, Loss1_betai2, Loss1_betai3, Loss1_betai4, Loss1_betai5,      Z      ; 
     Loss2_betai1, Loss2_betai2,      Z      , Loss2_betai4,      Z      , Loss2_betai6;
     Loss3_betai1, Loss3_betai2, Loss3_betai3, Loss3_betai4, Loss3_betai5,     Z       ;
     Loss4_betai1,      Z      , Loss4_betai3, Loss4_betai4, Loss4_betai5,     Z       ;
     Loss5_betai1, Loss5_betai2, Loss5_betai3, Loss5_betai4, Loss5_betai5, Loss5_betai6;
          Z      , Loss6_betai2,      Z      ,      Z      ,      Z      , Loss6_betai6];

end

