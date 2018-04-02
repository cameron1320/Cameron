function [ Kb,Cb ] = Bear(obj,kx,ky,c,DOF)
%% Creates bearing nodal stiffness matrix in the form:
% disp1y disp1z ang1y ang1z disp1x ang1x
% ______________________________________
%|                                      |disp1y
%|                                      |disp1z
%|                                      |ang1y
%|                                      |ang1z
%|                                      |disp1x
%|______________________________________|ang1x
%% Begin function
switch DOF
    case 6
        Kb = diag([kx,ky,0,0,0,0]);
        Cb = diag([c*kx,c*ky,0,0,0,0]);
    case 4
        Kb = diag([kx,ky,0,0]);
        Cb = diag([c*kx,c*ky,0,0]);
    case 2
        k = mean([kx,ky]);
        Kb = diag([k,k]);
        Cb = diag([c*k,c*k]);
end
end