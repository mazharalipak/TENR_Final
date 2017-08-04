%% Slack Variables.........

[ tnr ] = tnr_init(case9);



[pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);
 x = V2x_pf(tnr.V0, tnr);


V = x2V_pf(x, tnr);

Vmag_pq=V(pq);            % Complex voltages at PQ buses........

Va_pq=abs(Vmag_pq);            % Voltage magnitudes at PQ buses..............

Vm_pq=angle(Vmag_pq);          % Voltage angles at PQ buses...................


[pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);


%% Functions for V_min and V_max Feasibility Constarints..........




%%  partial derivative of V_min feasibility w.r.r v_a (angles)


J31=sparse(npq,npv+npq);

%% partial derivative of V_max feasibility w.r.r V_m (magnitudes)

J32=2*abs(V).*((cosd(angle(V))).^2+(sind(angle(V)).^2));
J32=sparse(diag(J32));
J32=J32(pq,pq);


%%

J41=J31;
J42=J32;



%% Slack Variables..............

Slack_H=ones(npq,1);

Slack_L=ones(npq,1);



%% partial derivative of V_min feasibility w.r.r Slack_L
 
J33=-2*Slack_L;                  % partial derivative of V_min feasibility w.r.r Slack_L
J33=sparse(diag(J33));

J34=sparse(npq,npq);             % partial derivative of V_min feasibility w.r.r Slack_H

%% partial derivative of V_max feasibility w.r.r Slack_H

J43=J34;                  % partial derivative of V_max feasibility w.r.r Slack_L

J44=2*Slack_H;            % partial derivative of V_max feasibility w.r.r Slack_H
J44=sparse(diag(J44));


%% Functions for V_min and V_max Feasibility Constarints..........

%Delta_V_max=abs(V(pq)) -(V_max*(ones(npq,1)).^2-(Slack_H).^2);

%Delta_V_min=abs(V(pq)) -(V_min*(ones(npq,1)).^2+(Slack_L).^2);


%% We need complex Volatge magnitudes........ and angles............













