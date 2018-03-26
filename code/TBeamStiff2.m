function [ K ] = TBeamStiff2( obj,E,I,l,r,nuv,poi,DOF )
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
I = A*r^2/4;
alph = 12*E*I/k/G/A/l^2;
syms x
z=x/l;

%%Shape Functions%%
N1 = 1-z;
N2 = z;
Tt1 = (1/(1+alph))*(2*z^3 - 3*z^2 - alph*z + 1 + alph);
Tt2 = (1/(1+alph))*(-2*z^3+3*z^2+alph*z);
Tr1 = (l/(1+alph))*(z^3 - (2+1/2*alph)*z^2 + (1+1/2*alph)*z);
Tr2 = (l/(1+alph))*(z^3 - (1-1/2*alph)*z^2 - 1/2*alph*z);
Rt1 = 6/l*(1/(1+alph))*(z^2-z);
Rt2 = -Rt1;
Rr1 = (1/(1+alph))*(3*z^2 - (4+alph)*z + 1 + alph);
Rr2 = (1/(1+alph))*(3*z^2 - (2-alph)*z);
%%%%%%%%%
t=DOF;
switch t
    case 6
        %Transformation Matricies%%
        D = [12*E*I/alph/l^2       0           0    0    0       0;
                    0         12*E*I/alph/l^2  0    0    0       0;
                    0               0           E*I 0   0       0;
                    0               0           0   E*I 0       0;
                    0               0           0   0   E*A     0;
                    0               0           0   0   0   12*E*I/alph/l^2/A*2*I];
        P = [0  0 0 0 0 0;
            0  0 0 0 0 0;
            0  1 0 0 0 0;
            -1  0 0 0 0 0;
            0  0 0 0 0 0;
            0  0 0 0 0 0];
        T = kron(diag([1,1,0]),[0,1;-1,0]);
        N = [Tt1  0    0  Tr1  Tt2  0    0  Tr2  0  0  0  0;
              0  Tt1  Tr1  0    0  Tt2  Tr2  0   0  0  0  0;
              0  Rt1  Rr1  0    0  Rt2  Rr2  0   0  0  0  0;
             Rt1  0    0  Rr1  Rt2  0    0  Rr2  0  0  0  0;
              0   0    0   0    0   0    0   0  N1  0 N2  0;
              0   0    0   0    0   0    0   0   0 N1  0 N2];
        I=diag([1,1,-1,1,1,1]);
        N = I*N*[I,zeros(length(I));zeros(length(I)),I];

    case 4
        D = [12*E*I/alph/l^2       0           0   0;
                  0           12*E*I/alph/l^2  0   0;
                  0                 0           E*I 0;
                  0                 0           0   E*I];
        P = [0 0 0 0;
            0 0 0 0;
             0 1 0  0;
             -1 0 0 0];
        T = kron(diag([1,1]),[0,1;-1,0]);
        N = [Tt1  0    0  Tr1  Tt2  0    0  Tr2;
              0  Tt1  Tr1  0    0  Tt2  Tr2  0;
              0  Rt1  Rr1  0    0  Rt2  Rr2  0;
             Rt1  0    0  Rr1  Rt2  0    0  Rr2];
        I=diag([1,1,-1,1]);  
        N = I*N*[I,zeros(length(I));zeros(length(I)),I];
    case 2
        D = [12*E*I/alph/l^2   0;
                    0           E*I];
        P = [0 0;
            -1 0];
        T =[0,1;-1,0];
        N = [Tt1  Tr1  Tt2  Tr2;
             Rt1  Rr1  Rt2  Rr2];
end





Bt = diff(N',x) + N'*P;
B=Bt';
K_B = int(B.'*D*B,x,0,l);
K_C = int(B.'*T*D*B,x,0,l);
K = K_B + nuv.*1i.*K_C;
if isa(l,'sym') || isa(E,'sym') || isa(I,'sym') || isa(nuv,'sym')
else
    K = double(K);
end
end

