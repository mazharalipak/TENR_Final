%% Normal Newton Raphson with Voltage Feasibility........
clc;
clear all;
step=1;
testcase = case9;
%% Loading Data structure for computing solutions..............

mpc = loadcase(testcase);
mpc = ext2int(mpc);
sol = runpf(mpc);
sol=ext2int(sol);

%% Input Data..........................

tnr = tnr_init(sol);
type=sol.bus(:,2);
x = V2x_pf(tnr.V0, tnr);

[pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);

V_max=1.1;
V_min=0.9;
Slack_H=ones(npq,1);
Slack_L=0.5*ones(npq,1);
l=0;
z = [];

Tole = 1;  
Iter = 1;
counter = 0;
while (Tole > 1e-3) 
    
[nx, ~] = deal(length(x), length(z));  
F = tnr.F(x, l, tnr);
[J,JJ] = tnr.J(x, l, tnr,sol);

%% Creating Jacobain with Voltage Feasibility Constraints...............

[Delta_V_max, Delta_V_min]=slack_delta(x, V_max,V_min,Slack_H,Slack_L,tnr);

[J31,J32,J33,J34,J41,J42,J43,J44]=J_pf_slack(x,Slack_L,Slack_H,tnr);

[aJ,bJ]=size(J);

J_maz=[J sparse(aJ,2*npq);J31 J32 J33 J34;J41 J42 J43 J44];
    
dMM=[F;Delta_V_min;Delta_V_max];
        
dv = - J_maz \ dMM;
    
dx = dv(1:nx);
dsl= dv(nx+1:nx+npq);
dsh= dv(nx+npq+1:end);

x  = x+step*dx;
dsl= Slack_L+step*dsl;
dsh= Slack_H+step*dsh;


      Iter = Iter + 1;
      Tole=max(abs(dMM));  
    counter = counter + 1;
    if counter ==50;
        break;
    end 
end

Tole
Iter