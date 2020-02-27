function H_on=H_onBDG(tw,uw,Delta,Lyw,Alphay,Vz)
tauz=[1,0;0,-1];taux=[0,1;1,0];tauy=[0,-1i;1i,0];
sgz=tauz;
H_on=kron(tauz,(4*tw-uw)*kron(eye(Lyw),eye(2)))+kron(eye(2),Vz*kron(eye(Lyw),sgz))+...
   kron(tauz, kron(diag(ones(Lyw-1,1),1),-tw*eye(2)+Alphay*(-1i/2)*sgz))+...
    kron(tauz,kron(diag(ones(Lyw-1,1),-1),-tw*eye(2)+Alphay*(1i/2)*sgz))+ kron(taux,kron(eye(Lyw),(Delta)*eye(2)));
end