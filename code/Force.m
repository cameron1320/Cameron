function F = Force(obj,md, Id, Ip, a, Ki,DOF)
%% Creates Disk nodal stiffness matrix in the form:
% 
% _ 
%| |disp1x
%| |disp1y
%| |ang1x
%|_|ang1y

%% Begin Function

Fx = md*a;
Fy = md*a;
Mx = (Id-Ip)*Ki;
My = (Id-Ip)*Ki;
% Fx = md*a*(w^2*cos(theta+phi1)-alph*sin(theta+phi1)); %N, Forcing funtion at disk 1
% Fy = md*a*(w^2*sin(theta+phi1)+alph*cos(theta+phi1)); %N, Forcing funtion at disk 2
% Mx = -w^2*(Ip - Id)*Ki*cos(theta+phi2); %N, Forcing funtion at disk 1
% My = w^2*(Ip - Id)*Ki*sin(theta+phi2); %N, Forcing funtion at disk 2;
switch DOF
    case 6
        F = [Fx; Fy; Mx; My; 0; 0];
    case 4
        F = [Fx; Fy; Mx; My];
    case 2
        F = [Fx; My];
end

end