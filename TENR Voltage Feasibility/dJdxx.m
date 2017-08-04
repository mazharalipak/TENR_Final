function [dJdxxx]=dJdxx(u,tnr)

[pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);

ux=u(1:npv+2*npq);

us1=u(npv+2*npq+1:npv+3*npq);
us2=u(npv+3*npq+1:end);
  
G32=sparse(diag(2*ux));
G32=G32(pq,pq);
G42=G32;

G33=sparse(diag(-2*us1));
G44=sparse(diag(2*us2));

dJdxxx=real([sparse(npq,npv+npq) G32 G33 sparse(npq,npq);sparse(npq,npv+npq) G42 sparse(npq,npq) G44]);
 
end