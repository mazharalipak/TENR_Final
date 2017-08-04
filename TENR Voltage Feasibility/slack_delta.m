function [Delta_V_max, Delta_V_min]=slack_delta(x, V_max,V_min,Slack_H,Slack_L,tnr)

V = x2V_pf(x, tnr);

[~, pq, ~, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);

Delta_V_max=((real(V(pq))).^2 + (imag(V(pq))).^2)-((V_max*(ones(npq,1))).^2 -(Slack_H).^2);

Delta_V_min=((real(V(pq))).^2 + (imag(V(pq))).^2)-((V_min*(ones(npq,1))).^2 +(Slack_L).^2);

end