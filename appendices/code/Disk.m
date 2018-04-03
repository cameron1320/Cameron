function [Md,Gd,Fd] = Disk(obj,md,Id,Ip,a,Ki,DOF)
%% Creates Disk nodal stiffness matrix in the form:
% disp1x disp1y ang1x ang1y
% _________________________ 
%|                         |disp1x
%|                         |disp1y
%|                         |ang1x
%|_________________________|ang1y

%% Begin Function
Fx = md*a;
Fy = md*a;
Mx = (Id-Ip)*Ki;
My = (Id-Ip)*Ki;
Gi=[0,Ip;-Ip,0];
switch DOF
    case 6
        Md = diag([md,md,Id,Id,md,Ip]);
        Gd = 1i*kron(diag([0,1,0]),Gi);
        Fd = [Fx; Fy; Mx; My; 0; 0];
    case 4
        Md = diag([md,md,Id,Id]);
        Gd = 1i*kron(diag([0,1]),Gi);
        Fd = [Fx; Fy; Mx; My];
    case 2
        Md = diag([md,Id]);
        Gd = 1i*[0,0;0,Ip];
        Fd = [Fx; My];
end
% Fx = md*a*(w^2*cos(theta+phi1)-alph*sin(theta+phi1)); %N, Forcing funtion at disk 1
% Fy = md*a*(w^2*sin(theta+phi1)+alph*cos(theta+phi1)); %N, Forcing funtion at disk 2
% Mx = -w^2*(Ip - Id)*Ki*cos(theta+phi2); %N, Forcing funtion at disk 1
% My = w^2*(Ip - Id)*Ki*sin(theta+phi2); %N, Forcing funtion at disk 2;
end