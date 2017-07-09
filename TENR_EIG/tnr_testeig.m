
function [Tole,l,Iter]=tnr_testeig(datain,step)
%% Input Test case Structure and Intitial start.........................
tic;
testcase = datain;
mpc = loadcase(testcase);
mpc = ext2int(mpc);
sol = runpf(mpc);
sol=ext2int(sol);

%% Input Data..........................

tnr = tnr_init(sol);
[~, ~, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);

x = V2x_pf(tnr.V0, tnr);
l = 1;
t=1;
z = 0.5*(ones((npv+2*npq),1));

%% Newtpn Iterations ................................
Tole = 1;  
Iter = 1;
counter = 0;

while (Tole > 1e-2)  
[ dx, dl,dz,dt,dMM] = tnr_step( x, l, z,t, tnr);
x=x+step*dx;
l=l+step*dl;
z=z+step*dz;
t=t+step*dt;

 Iter = Iter + 1;
      Tole=max(abs([dMM;dx;dl]));  
    counter = counter + 1;
    if counter ==100
        break;
    end    
end
toc;

end


