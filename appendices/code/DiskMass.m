function Md = DiskMass(obj,md,Id,ip,DOF)
%% Creates Disk nodal stiffness matrix in the form:
% disp1x disp1y ang1x ang1y
% _________________________ 
%|                         |disp1x
%|                         |disp1y
%|                         |ang1x
%|_________________________|ang1y

%% Begin Function
switch DOF
    case 6
        Md = diag([md,md,Id,Id,md,ip]);
    case 4
        Md = diag([md,md,Id,Id]);
    case 2
        Md = diag([md,Id]);
end

end