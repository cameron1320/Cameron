function [ A ] = Area( obj,do,di )
%AREA Summary of this function goes here
%   Detailed explanation goes here
if nargin == 4
    di = 0;
end

    A = pi*((do-di)/2)^2;        %m^2 cross sectional Area

end

