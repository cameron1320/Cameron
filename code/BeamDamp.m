function Ce = BeamDamp(obj,E,I,l,rho,r,Ip,nuv,poi)
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

syms x;
psi = [1 - 3*(x/l)^2 + 2*(x/l)^3, x - 2*l*(x/l)^2 + l*(x/l)^3, 3*(x/l)^2 - 2*(x/l)^3, -l*(x/l)^2 + l*(x/l)^3];
dpsi = diff(psi);
cpsi = [psi(1),      0,       0, psi(2), psi(3),      0,       0, psi(4);
            0, psi(1), -psi(2),      0,      0, psi(3), -psi(4),      0];
cphi = [      0, -dpsi(1), dpsi(2),       0,       0, -dpsi(3), dpsi(4),       0;
     dpsi(1),        0,       0, dpsi(2), dpsi(3),        0,       0, dpsi(4)];
shape = diff(diff(cpsi));
A = pi*r^2;
ip = 1/2*rho*A*r^2;
p = 24*(1+poi)*I/0.9/(pi*r^2)/l^2;
K_B = E*I*int(shape'*shape,x,0,l);
T = zeros(8);
T(3,3)=1;T(3,7)=-1;T(4,4)=1;T(4,8)=-1;T(7,3)=-1;T(7,7)=1;T(8,4)=-1;T(8,8)=1;
% K_B = (K_B + E*I/l*p.*T)./(1+p);
G = ip*int((cphi'*[0,1;-1,0]*cphi),x,0,l);
if isa(l,'sym') || isa(E,'sym') || isa(I,'sym') || isa(nuv,'sym')
else
K_B = double(K_B);
G = double(G);
end
Ce = nuv.*K_B - 1i.*G;

end

