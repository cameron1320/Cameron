function [ K_new ] = SCond(obj, K, c )
%SCOND uses static condensation to reduce the shear angle DOFs out of
%elemental matricies
%   Detailed explanation goes here
DOF = length(K);
for i = 1:1:DOF
    for j=1:1:DOF
        if i < c+1 && j < c+1
            k11(i,j) = K(i,j);
        end
        if i > c && j < c+1
            k12(i-c,j) = K(i,j);
        end
        if i > c && j > c
            k22(i-c,j-c) = K(i,j);
        end
    end
end
K_new = k11 - k12'*k22^-1*k12;         

end

