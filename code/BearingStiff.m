function Kb = BearingStiff(obj,kx,ky,DOF)
%% Creates bearing nodal stiffness matrix in the form:
% disp1x disp1y ang1x ang1y
% _________________________
%|                         |disp1x
%|                         |disp1y
%|                         |ang1x
%|_________________________|ang1y
%% Begin function
switch DOF
    case 6
        Kb = diag([kx,ky,0,0,0,0]);

    case 4
        Kb = diag([kx,ky,0,0]);

    case 2
        k = mean([kx,ky]);
        Kb = diag([k,k]);

end
end