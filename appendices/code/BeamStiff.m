function Ke = BeamStiff(obj,E,I,l,r,nuv,poi)
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
%
%   Update: 11-17-17
%           Include Shear correction term, for Timoshenko beam. p =
%           24*(1+poisson)*I/k_s/A/l^2, k_s=0.9 for solid circular shaft.

%% Begin Function
syms x;
psi = [1 - 3*(x/l)^2 + 2*(x/l)^3, x - 2*l*(x/l)^2 + l*(x/l)^3, 3*(x/l)^2 - 2*(x/l)^3, -l*(x/l)^2 + l*(x/l)^3];
cpsi = [psi(1),      0,       0, psi(2), psi(3),      0,       0, psi(4);
            0, psi(1), -psi(2),      0,      0, psi(3), -psi(4),      0];
shape = diff(diff(cpsi));
p = 24*(1+poi)*I/0.9/(pi*r^2)/l^2;
% syms p
K_B = E*I*int(shape'*shape,x,0,l);
T = zeros(8);
T(3,3)=1;T(3,7)=-1;T(4,4)=1;T(4,8)=-1;T(7,3)=-1;T(7,7)=1;T(8,4)=-1;T(8,8)=1;
% K_B = (K_B + E*I/l*p.*T)./(1+p);
% K_C = K_B*kron(eye(4),[0,1;-1,0]);
K_C = E*I*int(shape'*[0, 1;-1, 0]*shape,x,0,l);
if isa(l,'sym') || isa(E,'sym') || isa(I,'sym') || isa(nuv,'sym')
else
K_B = double(K_B);
K_C = double(K_C);
end
Ke = K_B + nuv.*1i.*K_C;


end