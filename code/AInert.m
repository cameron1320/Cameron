function [ I ] = AInert(obj, do, di )
%AINERT Summary of this function goes here
%   Detailed explanation goes here
if nargin == 4
    di = 0;
end

I = pi*(((do-di)/2)^4)/4;    %m^4 area moment of inertia

end

