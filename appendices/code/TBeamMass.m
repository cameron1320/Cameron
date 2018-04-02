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
Phi = 24*(1+poi)*I*(10/9)/(pi*r^2)/l^2;
syms zeta %A rho I
% syms l real
% zeta = z/l;
% syms Phi real
A = pi*r^2;

N11 = (1+Phi*(1-zeta)-3*zeta^2+2*zeta^3)/(1+Phi);
N12 = l*zeta*(1+(1/2)*Phi*(1-zeta)-2*zeta+zeta^2)/(1+Phi);
N13 = zeta*(Phi+3*zeta-2*zeta^2)/(1+Phi);
N14 = l*zeta*((-1/2)*Phi*(1-zeta)-zeta+zeta^2)/(1+Phi);
N21 = 6*zeta*(zeta-1)/(l*(1+Phi));
N22 = (1+Phi*(1-zeta)-4*zeta+3*zeta^2)/(1+Phi);
N23 = -6*zeta*(zeta-1)/(l*(1+Phi));
N24 = (Phi*zeta-2*zeta+3*zeta^2)/(1+Phi);

Nb = [N11   0    0 N12 N13   0    0 N14;
        0 N11 -N12   0   0 N13 -N14   0];
Ns = [N21   0    0 N22 N23   0    0 N24
        0 N21 -N22   0   0 N23 -N24   0];
% N3 = Ns - 1/l*diff(Nb,zeta);
% Nsd = diff(Ns,zeta);
% id = 1/4*rho*A*(r^2+x^2);
% id = 1/12*rho*A*l^3+1/4*rho*A*l*r^2;
% id=rho*I;
% ip=2*rho*I;
N = [Nb;Ns];
D = diag([rho*A*l,rho*A*l,Id,Id]);
Me = int(N'*D*N,zeta,0,1);

% Mt = rho*A*l.*int(Nb'*Nb,zeta,0,1);
% Mr = rho*I*l.*int(Ns'*Ns,zeta,0,1);

% Me = Mt + Mr;
if isa(rho,'sym') || isa(l,'sym') || isa(r,'sym')
else
Me = double(Me);
end


end