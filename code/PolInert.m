function [ Ip ] = PolInert( obj, do,l,rho,di )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 4
    di = 0;
end
 Ip = (1/2)*pi*rho*l*((do/2)-(di/2))^2*((do/2)^2 + (di/2)^2);

end

