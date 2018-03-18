function Md = DiskMass(obj,md,Id)
%% Creates Disk nodal stiffness matrix in the form:
% disp1x disp1y ang1x ang1y
% _________________________ 
%|                         |disp1x
%|                         |disp1y
%|                         |ang1x
%|_________________________|ang1y

%% Begin Function
Md = diag([md,md,Id,Id]);

end