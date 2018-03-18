function Kb = BearingStiff(obj,kx,ky)
%% Creates bearing nodal stiffness matrix in the form:
% disp1x disp1y ang1x ang1y
% _________________________
%|                         |disp1x
%|                         |disp1y
%|                         |ang1x
%|_________________________|ang1y
%% Begin function
Kb = diag([kx,ky,0,0]);
end