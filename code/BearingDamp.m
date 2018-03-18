function Cb = BearingDamp(obj,cx,cy)
%% Creates bearing nodal stiffness matrix in the form:
% disp1x disp1y ang1x ang1y
% _________________________
%|                         |disp1x
%|                         |disp1y
%|                         |ang1x
%|_________________________|ang1y
%% Begin function
Cb = diag([cx,cy,0,0]);

end