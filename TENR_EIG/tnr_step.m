function [ dx, dl,dz,dt,dMM] = tnr_step( x, l, z,t, tnr)
%tnr_step Summary of this function goes here
%   Detailed explanation goes here

%%  Routines for Creating Elements of Extended Jacobain.....

    [nx, ~] = deal(length(x), length(z));
    
    F = tnr.F(x, l, tnr);
    [J,JJ] = tnr.J(x, l, tnr);

    [G, u, v] = tnr.G(J, z);
    
    dGdz = tnr.dGdz(J, z);
    
    dGdxl = tnr.dJdxl(x, l, u, v, tnr);        
    %% enforcing Lambda>0 condition.............
    
    Delta_M=l-t^2;                                     % nonlinear equality equation for Lambda >0 ................
        
    ra=zeros(1,size(JJ,2)+1);
    ra(1,end-1)=1;
    ra(1,end)=-2*t;

%% Extended Jacobian Matrix................
        
    fullJ=[JJ zeros(size(JJ,1),1); dGdxl dGdz zeros(size(J,1)+1,1) zeros(size(dGdxl,1),1);ra];
        
      
    dMM=[F;G;Delta_M];
        
    dv = - fullJ \ [F; G;Delta_M];
    
    dx = dv(1:nx);
    dz = dv(nx+1:end-2);
    dl = dv(end-1);
    dt=  dv(end);
end

