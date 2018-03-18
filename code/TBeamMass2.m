function [Me] = TBeamMass(obj,rho,l,r,I,Id,Ip,poi)
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
% k = 6*(1+poi)/(7+6*poi);
% G = E/2/(1+poi);
A = pi*r^2;
% alph = 12*E*I/k/G/A/l^2;
alph = 24*(1+poi)*I*(10/9)/(pi*r^2)/l^2;

syms x
z=x/l;

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

% D = diag([rho*A,rho*A,rho*I,rho*I,rho*A,2*rho*I]);
D = diag([rho*A,rho*A,rho*I,rho*I]);

Me = int(N.'*D*N,x,0,l);

if isa(rho,'sym') || isa(l,'sym') || isa(r,'sym')
else
Me = double(Me);
end
% Me(9,:)=[];
% Me(9,:)=[];
% Me(9,:)=[];
% Me(9,:)=[];
% Me(:,9)=[];
% Me(:,9)=[];
% Me(:,9)=[];
% Me(:,9)=[];


end