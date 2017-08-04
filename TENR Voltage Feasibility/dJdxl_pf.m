function [ dJ] = dJdxl_pf(x, l, u, v,E_M, tnr)
%pf_dJ evaluate u'*dJ*v: sensitivity of the Jacobian
%   Detailed explanation goes here
    
    [pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);
    
    V = x2V_pf(x, tnr);
    
    ux=u(1:npv+npq+npq);
    
    lam = conj(F2S_pf(ux, tnr));
    
    [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(tnr.Ybus, V, lam);

%     dVa = v(1:npv+npq);
%     dVm_pq = v(npv+npq+1:npv+2*npq);
    
    H_dJ= real([Gaa([pv;pq],[pv;pq]) Gav([pv;pq],pq); 
        Gva(pq,[pv;pq])  Gvv(pq,pq)]);
    
    dJdxxx=dJdxx(u,tnr);
    
    [cj,~]=size(H_dJ);
    
    dJJ=[H_dJ sparse(cj,2*npq);dJdxxx]*E_M;
    dJ=dJJ'*v;
        
%     dJ = real([
%         (Gaa([pv;pq],[pv;pq]))*dVa + (Gav([pv;pq],pq))*dVm_pq; 
%         (Gva(pq,[pv;pq]))*dVa + (Gvv(pq,pq))*dVm_pq;
%         0
%     ]');

end

