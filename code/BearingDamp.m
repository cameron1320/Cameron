function Cb = BearingDamp(obj,cx,cy,DOF)
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
        Cb = diag([cx,cy,0,0,0,0]);
    case 4
        Cb = diag([cx,cy,0,0]);
    case 2
        c = mean([cx,cy]);
        Cb = diag([c,c]);
end


end