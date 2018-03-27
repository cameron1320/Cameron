function [ m ] = Mass(obj, do,l,rho,di )
%MASS Summary of this function goes here
%   Detailed explanation goes here
m = pi*((do-di)/2)^2*l*rho;

end

