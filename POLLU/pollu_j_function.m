function J = pollu_j_function(h_0, hd, ci, yi)
%% Chemical Parameters definition
% rate constants
k1 = 3.5e-1;
k2 = 2.66e1;
k3 = 1.23e4;
k4 = 8.6e-4;
k5 = 8.2e-4;
k6 = 1.5e4;
k7 = 1.3e-4;
k8 = 2.4e4;
k9 = 1.65e4;
k10 = 9e3;
k11 = 2.2e-2;
k12 = 1.2e4;
k13 = 1.88;
k14 =  1.63e4;
k15 = 4.8e6;
k16 = 3.5e-4;
k17 = 1.75e-2;
k18 = 1e8;
k19 = 4.44e11;
k20 = 1.24e3;
k21 = 2.1;
k22 = 5.78;
k23 = 4.74e-2;
k24 = 1.78e3;
k25 = 3.12;

Loss1_betai1 = ci*hd + (k1 + k10*yi(:, 11) + k14*yi(:, 6) + k23*yi(:, 4) + k24*yi(:, 19)).*h_0;
Loss1_betai2 = (- k2*yi(:, 4) - k3*yi(:, 5) - k9*yi(:, 11) - k12*yi(:, 10)).*h_0;
Z = zeros(size(h_0));
Loss1_betai4 = (k23*yi(:, 1) - k2*yi(:, 2)).*h_0;
Loss1_betai5 =  - k3*yi(:, 2).*h_0;
Loss1_betai6 = k14*yi(:, 1).*h_0;
Loss1_betai10 = - k12*yi(:, 2).*h_0;
Loss1_betai11 = (k10*yi(:, 1) - k9*yi(:, 2) ).*h_0;
Loss1_betai13 = - k11.*h_0;
Loss1_betai19 = (k24*yi(:, 1) - k22).*h_0;
Loss1_betai20 = - k25.*h_0;
        
Loss2_betai1 = - k1.*h_0;
Loss2_betai2 = ci*hd + (k2*yi(:, 4) + k3*yi(:, 5) + k9*yi(:, 11) + k12*yi(:, 10)).*h_0;
Loss2_betai4 = k2*yi(:, 2).*h_0;
Loss2_betai5 = k3*yi(:, 2).*h_0;
Loss2_betai10 = k12*yi(:, 2).*h_0;
Loss2_betai11 = k9*yi(:, 2).*h_0;
Loss2_betai19 = - k21.*h_0;

Loss3_betai1 = - k1.*h_0;
Loss3_betai3 = ci*hd + k15.*h_0;
Loss3_betai4 = - k17.*h_0;
Loss3_betai16 = - k19.*h_0;
Loss3_betai19 = - k22.*h_0;

Loss4_betai1 = k23.*h_0;
Loss4_betai2 = k2*yi(:, 4).*h_0;
Loss4_betai3 = - k15.*h_0;
Loss4_betai4 = ci*hd + (k2*yi(:, 2) + k16 + k17 + k23*yi(:, 1)).*h_0;

Loss5_betai2 = k3*yi(:, 5).*h_0;
Loss5_betai5 = ci*hd + k3*yi(:, 2).*h_0;
Loss5_betai6 = (- k6*yi(:, 7) - k20*yi(:, 17)).*h_0;
Loss5_betai7 = (- 2*k4 - k6*yi(:, 6)).*h_0;
Loss5_betai9 =  - k7.*h_0;
Loss5_betai14 = - k13.*h_0;
Loss5_betai17 = - k20*yi(:, 6).*h_0;

Loss6_betai1 = k14*yi(:, 6).*h_0;
Loss6_betai2 = - k3*yi(:, 5).*h_0;
Loss6_betai5 = - k3*yi(:, 2).*h_0;
Loss6_betai6 = ci*hd + (k6*yi(:, 7) + k8*yi(:, 9) + k14*yi(:, 1) + k20*yi(:, 17)).*h_0;
Loss6_betai7 = k6*yi(:, 6).*h_0;
Loss6_betai9 = k8*yi(:, 6).*h_0;
Loss6_betai16 = - 2*k18.*h_0;
Loss6_betai17 = k20.*h_0;

Loss7_betai6 = k6*yi(:, 7).*h_0;
Loss7_betai7 = ci*hd + (k4 + k5 + k6*yi(:, 6)).*h_0;
Loss7_betai14 = - k13.*h_0;

Loss8_betai6 = - k6*yi(:, 7).*h_0;
Loss8_betai7 = (- k4 - k5 - k6*yi(:, 6)).*h_0;
Loss8_betai8 = ci*hd;
Loss8_betai9 = - k7.*h_0;

Loss9_betai6 = k8*yi(:, 9).*h_0;
Loss9_betai9 = ci*hd + (k8*yi(:, 6) + k7).*h_0;

Loss10_betai2 = (- k9*yi(:, 11) + k12*yi(:, 10)).*h_0;
Loss10_betai9 = - k7.*h_0;
Loss10_betai10 = ci*hd + k12*yi(:, 2).*h_0;
Loss10_betai11 = - k9*yi(:, 2).*h_0;

Loss11_betai1 =  k10*yi(:, 11).*h_0;
Loss11_betai2 = k9*yi(:, 11).*h_0;
Loss11_betai6 = - k8*yi(:, 9).*h_0;
Loss11_betai9 = - k8*yi(:, 6).*h_0;
Loss11_betai11 = ci*hd + (k9*yi(:, 2) + k10*yi(:, 1)).*h_0;
Loss11_betai13 = - k11.*h_0;

Loss12_betai2 = - k9*yi(:, 11).*h_0;
Loss12_betai11 = - k9*yi(:, 2).*h_0;
Loss12_betai12 = ci*hd;

Loss13_betai1 =  - k10*yi(:, 11).*h_0;
Loss13_betai2 = - k10*yi(:, 1).*h_0;
Loss13_betai13 = ci*hd + k11.*h_0;

Loss14_betai2 = - k12*yi(:, 10).*h_0;
Loss14_betai10 = - k12*yi(:, 2).*h_0;
Loss14_betai14 = ci*hd + k13.*h_0;

Loss15_betai1 =  - k14*yi(:, 6).*h_0;
Loss15_betai6 = - k14*yi(:, 1).*h_0;
Loss15_betai15 = ci*hd;
 
Loss16_betai4 = - k16.*h_0;
Loss16_betai16 = ci*hd + (k18 + k19).*h_0;

Loss17_betai6 = k20*yi(:, 17).*h_0;
Loss17_betai17 = ci*hd + k20*yi(:, 6).*h_0;       

Loss18_betai6 = - k20*yi(:, 17).*h_0;
Loss18_betai17 = - k20*yi(:, 6).*h_0;
Loss18_betai18 = ci*hd;       

Loss19_betai1 = (k24*yi(:, 19) - k23*yi(:, 4)).*h_0;
Loss19_betai4 = - k23*yi(:, 1).*h_0;
Loss19_betai19 = ci*hd + (k21 + k22 + k24*yi(:, 1)).*h_0;
Loss19_betai20 = - k25.*h_0;

Loss20_betai1 =  - k24*yi(:, 19).*h_0;
Loss20_betai19 = - k24*yi(:, 1).*h_0;
Loss20_betai20 = ci*hd + k25.*h_0;


J = [Loss1_betai1  , Loss1_betai2  ,       Z       , Loss1_betai4  , Loss1_betai5  , Loss1_betai6  ,       Z       , ...
           Z       ,       Z       , Loss1_betai10 , Loss1_betai11 ,       Z       , Loss1_betai13 ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       , Loss1_betai19 , Loss1_betai20 ;
     Loss2_betai1  , Loss2_betai2  ,       Z       , Loss2_betai4  , Loss2_betai5  ,       Z       ,       Z       , ...
           Z       ,       Z       , Loss2_betai10 , Loss2_betai11 ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       , Loss2_betai19 ,       Z       ;
     Loss3_betai1  ,       Z       , Loss3_betai3  , Loss3_betai4  ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       , Loss3_betai16 ,       Z       ,       Z       , Loss3_betai19 ,       Z       ;
     Loss4_betai1  , Loss4_betai2  , Loss4_betai3  , Loss4_betai4  ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ;
           Z       , Loss5_betai2  ,       Z       ,       Z       , Loss5_betai5  , Loss5_betai6  , Loss5_betai7  , ...
           Z       , Loss5_betai9  ,       Z       ,       Z       ,       Z       ,       Z       , Loss5_betai14 , ...
           Z       ,       Z       , Loss5_betai17 ,       Z       ,       Z       ,       Z       ;
     Loss6_betai1  , Loss6_betai2  ,       Z       ,       Z       , Loss6_betai5  , Loss6_betai6  , Loss6_betai7  , ...
           Z       , Loss6_betai9  ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       , Loss6_betai16 , Loss6_betai17 ,       Z       ,       Z       ,       Z       ;
           Z       ,       Z       ,       Z       ,       Z       ,       Z       , Loss7_betai6  , Loss7_betai7  , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , Loss7_betai14 , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ;
           Z       ,       Z       ,       Z       ,       Z       ,       Z       , Loss8_betai6  , Loss8_betai7  , ...
     Loss8_betai8  , Loss8_betai9  ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ;
           Z       ,       Z       ,       Z       ,       Z       ,       Z       , Loss9_betai6  ,       Z       , ...
           Z       , Loss9_betai9  ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ;
           Z       , Loss10_betai2 ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       , Loss10_betai9 , Loss10_betai10, Loss10_betai11,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ;
     Loss11_betai1 , Loss11_betai2 ,       Z       ,       Z       ,       Z       , Loss11_betai6 ,       Z       , ...
           Z       , Loss11_betai9 ,       Z       , Loss11_betai11,       Z       , Loss11_betai13,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ;
           Z       , Loss12_betai2 ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       , Loss12_betai11, Loss12_betai12,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ;
     Loss13_betai1 , Loss13_betai2 ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       , Loss13_betai13,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ;
           Z       , Loss14_betai2 ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       , Loss14_betai10,       Z       ,       Z       ,       Z       , Loss14_betai14, ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ;
     Loss15_betai1 ,       Z       ,       Z       ,       Z       ,       Z       , Loss15_betai6 ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
     Loss15_betai15,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ;
           Z       ,       Z       ,       Z       , Loss16_betai4 ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       , Loss16_betai16,       Z       ,       Z       ,       Z       ,       Z       ;
           Z       ,       Z       ,       Z       ,       Z       ,       Z       , Loss17_betai6 ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       , Loss17_betai17,       Z       ,       Z       ,       Z       ;
           Z       ,       Z       ,       Z       ,       Z       ,       Z       , Loss18_betai6 ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       , Loss18_betai17, Loss18_betai18,       Z       ,       Z       ;
     Loss19_betai1 ,       Z       ,       Z       , Loss19_betai4 ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       , Loss19_betai19, Loss19_betai20;
     Loss20_betai1 ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       ,       Z       , ...
           Z       ,       Z       ,       Z       ,       Z       , Loss20_betai19, Loss20_betai20];

end

