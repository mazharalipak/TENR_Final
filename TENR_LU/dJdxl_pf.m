function [ dJ] = dJdxl_pf(x, l, u, v, tnr)
%pf_dJ evaluate u'*dJ*v: sensitivity of the Jacobian
%   Detailed explanation goes here
    
    [pv, pq, ~, ~]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);
    
    V = x2V_pf(x, tnr);
    
    
    lam = conj(F2S_pf(u, tnr));
    
    [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(tnr.Ybus, V, lam);
    
    H_dJ= [Gaa([pv;pq],[pv;pq]) Gav([pv;pq],pq); 
        Gva(pq,[pv;pq])  Gvv(pq,pq)];
        
    dJ= real([H_dJ*v;0]');
    
end

