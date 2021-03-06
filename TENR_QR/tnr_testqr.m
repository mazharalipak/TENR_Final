
function [Tole,Iter,l]=tnr_testqr(datain,step)
%% Input Test case Structure and Intitial start.........................
t=1;                              % Slack variable for \lambda >0 .....................

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
z = [];

%% Newtpn Iterations ................................
tic;

Tole = 1;  
Iter = 1;
counter = 0;
while (Tole > 1e-2)  
 
%% TENR Iteration 
   
%%
[ dx, dl,dt,dMM,~] = tnr_step( x, l,z,t, tnr);

%% Updating state variables......

%[~, G, ~] = svds(J, 1,0);

x=x+step*dx;
l=l+step*dl;
t=t+step*dt;

%% Storing stability index............

% stab_index(Iter)=G;
% rel_stab_index=stab_index./stab_index(1);
 
%% Checking Tolerance.......

 Iter = Iter + 1;
      Tole=max(abs([dMM;dx;dl]));  
    counter = counter + 1;
    if counter ==40
        break;
    end    
end


toc;

% semilogy(rel_stab_index)
% xlim([-Inf Inf])
% ylim([-Inf Inf])
end



