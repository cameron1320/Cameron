function [ K ] = TBeamStiff2( obj,E,I,l,r,nuv,poi )
%UNTITLED2 Summary of this function goes here
% disp1x disp1y ang1x ang1y disp2x disp2y ang2x ang2y
% ___________________________________________________
%|                                                   |disp1x
%|                                                   |disp1y
%|                                                   |ang1x
%|                                                   |ang1y
%|                                                   |disp2x
%|                                                   |disp2y
%|                                                   |ang2x
%|___________________________________________________|ang2y
k = 6*(1+poi)/(7+6*poi);
G = E/2/(1+poi);
A = pi*r^2;
alph = 12*E*I/k/G/A/l^2;
syms x
z=x/l;
%%Transformation Matricies%%
% P = [0  0 0 0 0 0;
%      0  0 0 0 0 0;
%      0 -1 0 0 0 0;
%      1  0 0 0 0 0;
%      0  0 0 0 0 0;
%      0  0 0 0 0 0];
% D = [k*G*A     0   0   0   0       0;
%          0 k*G*A   0   0   0       0;
%          0     0 E*I   0   0       0;
%          0     0   0 E*I   0       0;
%          0     0   0   0 E*A       0;
%          0     0   0   0   0 k*G*2*I];
% T = kron(diag([1,1,0]),[0,1;-1,0]);
 D = [12*E*I/alph/l^2    0   0   0;
         0   12*E*I/alph/l^2 0   0;
         0     0  E*I  0;
         0     0   0 E*I];
P = [0  0 0 0;
     0  0 0 0;
     0 -1 0 0;
     1  0 0 0];
T = kron(diag([1,1]),[0,1;-1,0]);

N1 = 1-z;
N2 = z;
Tt1 = (1/(1+alph))*(2*z^3 - 3*z^2 - alph*z + 1 + alph);
Tt2 = (1/(1+alph))*(-2*z^3+3*z^2-alph*z);
Tr1 = (l/(1+alph))*(z^3 - (2-1/2*alph)*z^2 + (1-1/2*alph)*z);
Tr2 = (l/(1+alph))*(z^3 - (1+1/2*alph)*z^2 + 1/2*alph*z);
Rt1 = 6/l*(1/(1+alph))*(z^2-z);
Rt2 = -Rt1;
Rr1 = (l/(1+alph))*(3*z^2 - (4-alph)*z + 1 - alph);
Rr2 = (l/(1+alph))*(3*z^2 - (2+alph)*z);

% N = [Tv1   0    0 Tt1 Tv2   0    0 Tt2  0  0  0  0;
%        0 Tv1  Tt1   0   0 Tv2  Tt2   0  0  0  0  0;
%        0   R  Rt1   0   0   -R  Rt2   0  0  0  0  0;
%        R   0    0 Rt1  -R   0    0 Rt2  0  0  0  0;
%        0   0    0   0   0   0    0   0 N1  0 N2  0;
%        0   0    0   0   0   0    0   0  0 N1  0 N2];
N = [Tt1   0    0 Tr1 Tt2   0    0 Tr2;
       0 Tt1  -Tr1   0   0 Tt2  -Tr2   0;
       0   Rt1  Rr1   0   0   Rt2  Rr2   0;
       -Rt1   0    0 Rr1  -Rt2   0    0 Rr2];
B = diff(N,x) - P.'*N;

K_B = int(B.'*D*B,x,0,l);
K_C = int(B.'*T*D*B,x,0,l);

K = K_B + nuv.*1i.*K_C;
if isa(l,'sym') || isa(E,'sym') || isa(I,'sym') || isa(nuv,'sym')
else
K = double(K);
end
% K(9,:)=[];
% K(9,:)=[];
% K(9,:)=[];
% K(9,:)=[];
% K(:,9)=[];
% K(:,9)=[];
% K(:,9)=[];
% K(:,9)=[];
end

