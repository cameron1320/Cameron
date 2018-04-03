function [Me] = TBeamMass(obj,rho,l,r,I,Id,Ip,poi,E,DOF)
%% Creates beam mass matrix in the form:
% disp1x disp1y ang1x ang1y disp2x disp2y ang2x ang2y
% ___________________________________________________
%|                                                   |disp1
%|                                                   |disp1y
%|                                                   |ang1x
%|                                                   |ang1y
%|                                                   |disp2x
%|                                                   |disp2y
%|                                                   |ang2x
%|___________________________________________________|ang2y
%% Begin Function
k = 6*(1+poi)/(7+6*poi);
G = E/2/(1+poi);
A = pi*r^2;
% I = A*r^2/4;
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
        Mi = diag([rho*A,rho*A,Id,Id,rho*A,2*rho*I]);
        N = [Tt1  0    0  Tr1  Tt2  0    0  Tr2  0  0  0  0;
              0  Tt1  Tr1  0    0  Tt2  Tr2  0   0  0  0  0;
              0  Rt1  Rr1  0    0  Rt2  Rr2  0   0  0  0  0;
             Rt1  0    0  Rr1  Rt2  0    0  Rr2  0  0  0  0;
              0   0    0   0    0   0    0   0  N1  0 N2  0;
              0   0    0   0    0   0    0   0   0 N1  0 N2];
         I=diag([1,1,-1,1,1,1]);
        N = I*N*[I,zeros(length(I));zeros(length(I)),I];
    case 4
        Mi = diag([rho*A,rho*A,rho*I,rho*I]);
        N = [Tt1  0    0  Tr1  Tt2  0    0  Tr2;
              0  Tt1  Tr1  0    0  Tt2  Tr2  0;
              0  Rt1  Rr1  0    0  Rt2  Rr2  0;
             Rt1  0    0  Rr1  Rt2  0    0  Rr2];
         I=diag([1,1,-1,1]);
        N = I*N*[I,zeros(length(I));zeros(length(I)),I];
    case 2
        Mi = diag([rho*A,rho*I]);
        N = [Tt1  Tr1  Tt2  Tr2;
             Rt1  Rr1  Rt2  Rr2];

end


Me = int(N.'*Mi*N,x,0,l);

if isa(rho,'sym') || isa(l,'sym') || isa(r,'sym')
else
Me = double(Me);
end


end