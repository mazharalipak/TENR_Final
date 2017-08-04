function [ dx,dsl,dsh,dl,dMM] = tnr_step( x, l,z, tnr,sol,V_max,V_min,Slack_H,Slack_L)
%tnr_step Summary of this function goes here
%   Detailed explanation goes here


%%  Routines for Creating Elements of Extended Jacobain.....

[nx, ~] = deal(length(x), length(z));
    
F = tnr.F(x, l, tnr);
[J,JJ] = tnr.J(x, l, tnr,sol);

[pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);
%% Creating Jacobain with Voltage Feasibility Constraints...............

[Delta_V_max, Delta_V_min]=slack_delta(x, V_max,V_min,Slack_H,Slack_L,tnr);
[J31,J32,J33,J34,J41,J42,J43,J44]=J_pf_slack(x,Slack_L,Slack_H,tnr);

[aJ,~]=size(J);

J_maz=[J sparse(aJ,2*npq);J31 J32 J33 J34;J41 J42 J43 J44];

[ G, u, v,E_M] = tnr.G(J_maz, z);
dGdxl = tnr.dJdxl(x, l, u, v,E_M, tnr);
    
        
%% Extended Jacobian Matrix................

aJ_maz=[JJ(:,end);zeros(2*npq,1)];

fullJJ=[J_maz aJ_maz;dGdxl' 0];
                             
%% Solving Linear system of Equation Newton Step....

dMM=[F;Delta_V_min;Delta_V_max;G];
        
dv = - fullJJ \ dMM;
    
dx = dv(1:nx);
dsl= dv(nx+1:nx+npq);
dsh= dv(nx+npq+1:end-1);
dl= dv(end);
end

