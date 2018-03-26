function Ce = TBeamDamp(obj,E,I,l,rho,r,Ip,nuv,poi)
%% Creates beam stiffness matrix in the form:
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

Phi = 24*(1+poi)*I*(10/9)/(pi*r^2)/l^2;
syms zeta
% syms l real
% zeta = z/l;
% syms Phi real
A = pi*r^2;
% ip = rho*A*l*r^2/2;

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
Nsd = diff(Ns,zeta);
N = [Nb;Ns];
DG = kron(diag([0,1]),[0, Ip;-Ip,0]);
N3 = Ns - 1/l*diff(Nb,zeta);
NK = [Nsd;N3];
DK = diag([E*I/l,E*I/l,12*E*I/l/Phi,12*E*I/l/Phi]);
% DK = diag([12*E*I/l/Phi,12*E*I/l/Phi,E*I/l,E*I/l]);

% Mr = rho*I*l.*int(Ns'*[0 1;-1 0]*Ns,zeta,0,1);
G = int(N'*DG*N,zeta,0,1);
K_B = int(NK'*DK*NK,zeta,0,1);


if isa(l,'sym') || isa(E,'sym') || isa(I,'sym') || isa(nuv,'sym')
else
K_B = double(K_B);
G = double(G);
end
Ce = nuv.*K_B - 1i.*G;

end

