function [ dJ,Gaa, Gav, Gva, Gvv] = dJdxl_pf( x, l, u, v, tnr)
%pf_dJ evaluate u'*dJ*v: sensitivity of the Jacobian
%   Detailed explanation goes here
    
    [pv, pq, ~, ~]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);
    
    V = x2V_pf(x, tnr);
        
    lam = conj(F2S_pf(v, tnr));
    
    [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(tnr.Ybus, V, lam);
    
    H_dJ= real([Gaa([pv;pq],[pv;pq]) Gav([pv;pq],pq); 
        Gva(pq,[pv;pq])  Gvv(pq,pq)]);    
   
    
    dJ=[H_dJ;zeros(1,size(H_dJ,1))];

end

