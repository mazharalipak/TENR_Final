function [ G, u, v,E_M] = G_svd(J_maz, z)
%tr_svd Calculate the value of G and its derivatives
%   x - variable 
%   J - Jacobian function: evaluates J(x)

 [u, G, v] = svds(J_maz, 1,0);
 E_M=sparse(eye(size(J_maz)));        %  I matrix

end

