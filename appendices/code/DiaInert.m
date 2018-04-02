function [ Id ] = DiaInert(obj, do,l,rho,di )
%DIAINERT Summary of this function goes here
%   Detailed explanation goes here
if nargin == 4
    di = 0;
end
Id = pi*rho*l/12*(3*((do/2)^4 - (di/2)^4) + l^2*((do/2)^4 - (di/2)^4));

end

