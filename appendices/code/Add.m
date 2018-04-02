function [Mr] = Add(obj, M, Me, index)
%% Add adds an element stiffness matrix to the global stiffness matrix
%       Mr == Return Matrix with the new element
%       M  == Existing global matrix before new element addition
%       Me == New element stiffness matrix
%       ei == Node locations corresponding to the new element stiffness
%       matrix

%% function %%

if isa(Me, 'sym')
    Mr = sym(M);
else
Mr = M;
end
m = length(index);
for i = 1:m
    for j = 1:m
        Mr(index(i),index(j)) = M(index(i),index(j)) + Me(i,j);
    end
end
end