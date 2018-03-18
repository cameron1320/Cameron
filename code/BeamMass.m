function [Me] = BeamMass(obj,rho,l,r,I,Id,Ip,poi)
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
syms x;
psi = [1 - 3*(x/l)^2 + 2*(x/l)^3, x - 2*l*(x/l)^2 + l*(x/l)^3, 3*(x/l)^2 - 2*(x/l)^3, -l*(x/l)^2 + l*(x/l)^3];
dpsi = diff(psi);
cpsi = [psi(1),      0,       0, psi(2), psi(3),      0,       0, psi(4);
            0, psi(1), -psi(2),      0,      0, psi(3), -psi(4),      0];
cphi = [      0, -dpsi(1), dpsi(2),       0,       0, -dpsi(3), dpsi(4),       0;
     dpsi(1),        0,       0, dpsi(2), dpsi(3),        0,       0, dpsi(4)];
A = pi*r^2;
% id = 1/4*rho*A*(r^2+x^2);
% id = 1/12*rho*A*l^3+1/4*rho*A*l*r^2;
% id=rho*I;
% ip=2*rho*I;
Mt = Id.*int(cpsi'*cpsi,x,0,l);
Mr = Ip.*int((cphi'*cphi),x,0,l);

Me = Mt + Mr;
if isa(rho,'sym') || isa(l,'sym') || isa(r,'sym')
else
Me = double(Me);
end


end