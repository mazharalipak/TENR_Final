function [J31,J32,J33,J34,J41,J42,J43,J44]=J_pf_slack(x,Slack_L,Slack_H,tnr)

[pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);

V = x2V_pf(x, tnr);
  
 %%  partial derivative of V_min feasibility w.r.r v_a
 
J31=sparse(npq,npv+npq);

%% partial derivative of V_max feasibility w.r.r V_m (magnitudes)

J32=2*abs(V);
J32=sparse(diag(J32));
J32=J32(pq,pq);

%%

J41=J31;
J42=J32;

%% partial derivative of V_min feasibility w.r.r Slack_L
 
J33=-2*Slack_L;                  % partial derivative of V_min feasibility w.r.r Slack_L
J33=sparse(diag(J33));

J34=sparse(npq,npq);             % partial derivative of V_min feasibility w.r.r Slack_H

%% partial derivative of V_max feasibility w.r.r Slack_H

J43=J34;                         % partial derivative of V_max feasibility w.r.r Slack_L

J44=2*Slack_H;                   % partial derivative of V_max feasibility w.r.r Slack_H
J44=sparse(diag(J44));


end