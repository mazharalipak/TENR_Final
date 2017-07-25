function [ G, u, v] = G_lu(J, z) 
%tr_lu Calculate the value of G and its derivatives
%   x - variable 
%   J - Jacobian function: evaluates J(x)

%% First option for this LU Transversality is based on a self written lu decomposition rather than Matlab

tol=1e-3;
[L_Mat, U_Mat, P_Mat, Q_Mat]=lucp(J,tol,'sparse');  

New_UMat=bsxfun(@rdivide, U_Mat(1:end,:), diag(U_Mat));

b=zeros(size(J,1),1);                            %% Right hand side ....(nth column of Identity matrix)
b(end,1)=1;

v=(New_UMat)\b;                                  %% The Right eigenvector in g_LU Transversality

u=((L_Mat)')\b;                                  %% The Left eigenvector in g_LU Transversality

G=U_Mat(end,end);                                %% g_LU
  
%% First option for this LU Transversality is based on Matlab LU decomposition.........................

% [L_Mat, U_Mat, P_Mat, Q_Mat]=lu(J); 
% New_UMat=bsxfun(@rdivide, U_Mat(1:end,:), diag(U_Mat));

% [aa,bb]=min(abs(diag(U_Mat)));
 
% b=zeros(size(J,1),1);                %% Right hand side ....(nth column of Identity matrix)
% b(bb,1)=1;
 
% v=(New_UMat)\b;                      %% The Right eigenvector in g_LU Transversality             
 
% u=((L_Mat))'\b;                      %% The Left eigenvector in g_LU Transversality  
 
% G=U_Mat(bb,bb);                      %% g_LU


end

