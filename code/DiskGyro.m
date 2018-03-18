function [G] = DiskGyro(obj,Ip)
%% Creates Disk nodal stiffness matrix in the form:
% disp1x disp1y disp1z ang1x ang1y
% ________________________________
%|                                |disp1x
%|                                |disp1y
%|                                |disp1z
%|                                |ang1x
%|________________________________|ang1y
%% Begin Function

syms w;
G = [0,0,0,0;0,0,0,0;0,0,0,1i*Ip;0,0,-1i*Ip,0];
end