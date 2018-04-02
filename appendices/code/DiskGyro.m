function [G] = DiskGyro(obj,Ip,DOF)
%% Creates Disk nodal stiffness matrix in the form:
% disp1x disp1y disp1z ang1x ang1y
% ________________________________
%|                                |disp1x
%|                                |disp1y
%|                                |disp1z
%|                                |ang1x
%|________________________________|ang1y
%% Begin Function
Gi=[0,-Ip;Ip,0];
switch DOF
    case 6
        G = -1i*kron(diag([0,1,0]),Gi);
    case 4
        G = -1i*kron(diag([0,1]),Gi);
    case 2
        G = 1i*[0,0;0,Ip];
end
end