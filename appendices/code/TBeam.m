function [ K,D,M ] = TBeam( obj, E,I,l,A,rho,Id,Ip,nuv,poi,DOF )
%TBeam Timoshenko beam element
% disp1y disp1z ang1y ang1z disp2y disp2z ang2y ang2z disp1x ang1x disp2x ang2x
% ____________________________________________________________________________
%|                                                                            |disp1y
%|                                                                            |disp1z
%|                                                                            |ang1y
%|                                                                            |ang1z
%|                                                                            |disp2y
%|                                                                            |disp2z
%|                                                                            |ang2y
%|                                                                            |ang2z
%|____________________________________________________________________________|ang2x
%% Properties %%
k = 6*(1+poi)/(7+6*poi);
G = E/2/(1+poi);
alph = 12*E*I/k/G/A/l^2;
%% Variables %%
syms x
z = x/l;
%% Shape Functions %%
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
%% Transformation Matricies %%
P = [0  0 0 0 0 0;
    0  0 0 0 0 0;
    0 -1 0 0 0 0;
    1  0 0 0 0 0;
    0  0 0 0 0 0;
    0  0 0 0 0 0];
Di = [12*E*I/alph/l^2          0       0   0   0           0;
    0          12*E*I/alph/l^2 0   0   0           0;
    0                  0       E*I 0   0           0;
    0                  0       0   E*I 0           0;
    0                  0       0   0   E*A         0;
    0                  0       0   0   0   12*E*I/alph/l^2/A*2*I];
Mi = diag([rho*A*l,rho*A*l,Id,Id,rho*A*l,Ip])./l;
Gi = kron(diag([0,1,0]),[0, Ip;-Ip,0])./l;
T = kron(diag([1,1,0]),[0,1;-1,0]);
N = [   Tt1   0    0  Tr1 Tt2  0    0  Tr2 0  0  0  0;
    0   Tt1  Tr1  0   0  Tt2  Tr2  0  0  0  0  0;
    0   Rt1  Rr1  0   0  Rt2  Rr2  0  0  0  0  0;
    Rt1   0    0  Rr1 Rt2  0    0  Rr2 0  0  0  0;
    0    0    0   0   0   0    0   0  N1 0  N2 0;
    0    0    0   0   0   0    0   0  0  N1 0  N2];
%% DOF Modifications%%
switch DOF
    case 6
        I = diag([1,1,-1,1,1,1]);
        N = I*N*[I,zeros(6);zeros(6),I]; % Adjust for use of -angy,-ang1,2y definitions in shape functions
    case 4
        P = P(1:4,1:4);
        Di = Di(1:4,1:4);
        Mi = Mi(1:4,1:4);
        Gi = Gi(1:4,1:4);
        T = T(1:4,1:4);
        I = diag([1,1,-1,1]); % Adjust for use of -angy,-ang1,2y definitions in shape functions
        N = I*N(1:4,1:8)*[I,zeros(4);zeros(4),I];
    case 2
        P = P(2:3,2:3); 
        Di = Di(2:3,2:3);
        Mi = Mi(2:3,2:3);
        Gi = Gi(2:3,2:3);
        T = T(2:3,2:3);
        N = N(2:3,[2,3,6,7]);
end
%% Elemental Integration %%
B = diff(N.',x) - N.'*P;
B = B.';
K_B = int(B.'*Di*B,x,0,l);
K_C = int(B.'*T*Di*B,x,0,l);
G = int(N.'*Gi*N,x,0,l);
M = int(N.'*Mi*N,x,0,l);
%% Case of symbolic Integration %%
if isa(l,'sym') || isa(E,'sym') || isa(I,'sym') || isa(nuv,'sym')
else
    K_B = double(K_B);
    K_C = double(K_C);
    G = double(G);
    M = double(M);
end
K =      K_B + nuv.*1i.*K_C;
D = nuv.*K_B +      1i.*G;

end