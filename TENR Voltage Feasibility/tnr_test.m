
function [Tole,l,Iter,x,dMM]=tnr_test(datain,step,va,vb)
%% Initial seed / germ for starting TENR Algorithm...........................

testcase = datain;

%% Loading Data structure for computing solutions..............

mpc = loadcase(testcase);
mpc = ext2int(mpc);
sol = runpf(mpc);
sol=ext2int(sol);

%% Input Data..........................

tnr = tnr_init(sol);
type=sol.bus(:,2);
sol.bus(:,2)=type;
x = V2x_pf(tnr.V0, tnr);
l = 1;

[~,~,~, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);

% V_max=abs(tnr.V0(1))+0.05*(abs(tnr.V0(1)));
% V_min=abs(tnr.V0(1))-0.05*(abs(tnr.V0(1)));

%% Voltage magnitudes constraints .............................

V_max=va;
V_min=vb;

Slack_H=ones(npq,1);
Slack_L=0.5*ones(npq,1);
z = [];

%% Newtpn Iterations ................................


tic;
Tole = 1;  
Iter = 1;
counter = 0;
while (Tole > 1e-2)  
   
%%.............

[ dx,dsl,dsh,dl,dMM] = tnr_step( x, l,z, tnr,sol,V_max,V_min,Slack_H,Slack_L);
%% Updating state variables......

x=x+step*dx;
Slack_L= Slack_L+step*dsl;
Slack_H= Slack_H+ step*dsh;
l=l+step*dl;

 
%% Checking Tolerance.......

 Iter = Iter + 1;
      Tole=max(abs(dMM));  
    counter = counter + 1;
    if counter ==150;
        break;
    end    
end
toc;
end



