function gw=Andreev_end(E,eta,tw,uw_Au,uw_EuS,V_x1,V_EuS,V_Au,Alphax,Alphay,Delta,Lx1,Lx2,Lx3, Lyw1,Lyw2,Lyw,L_y1,G_0n,tc_G,L_leadx,Gamma)
sgy=[0,-1i;1i,0];
sgx=[0,1;1,0];
sgz=[1,0;0,-1];
taux=sgx;
tauz=sgz;
tau0=eye(2);
i=1:Lx2;
Z0=Gamma/sqrt(Delta^2-(E+1i*eta)^2);
Delta0=Z0*Delta;
u=(1-i/(Lx2))*uw_Au+i/(Lx2)*uw_EuS;
V_u=(1-i/(Lx2))*V_Au+i/(Lx2)*V_EuS;
v=fliplr(u);
VV=flip(V_u);
%   HBDG=kron(tauz,kron(speye(Lxw),(4*tw)*kron(speye(Lyw),eye(2))+...
%         kron(diag(ones(Lyw-1,1),1),-tw*eye(2)+Alphay*(-1i/2)*sgz)+...
%         kron(diag(ones(Lyw-1,1),-1),-tw*eye(2)+Alphay*(1i/2)*sgz))+...
%         kron(diag(ones(Lxw-1,1),1),kron(speye(Lyw),-tw*eye(2)-Alphax*(-1i/2)*sgy))+...
%         kron(diag(ones(Lxw-1,1),-1),kron(speye(Lyw),-tw*eye(2)-Alphax*(1i/2)*sgy)))+...
%         Z*H_u+kron(tau0,H_z)+kron(taux,kron(speye(Lyw*Lxw),(Delta_r)*eye(2)));


Vc_BDG=kron(tauz,kron(eye(Lyw),-tw*eye(2)-Alphax*(-1i/2)*sgy));
VcL_B=kron(tauz,0*tc_G*kron(eye(Lyw),eye(2)));

GnnR=inv((E+1i*eta)*(1+Z0)*eye(4*Lyw)-H_onBDG(tw,uw_Au,Delta0,Lyw,Alphay,(1+Z0)*V_Au)-VcL_B*G_0n*VcL_B');
GnnL=inv((E+1i*eta)*(1+Z0)*eye(4*Lyw)-H_onBDG(tw,uw_Au,Delta0,Lyw,Alphay,(1+Z0)*V_Au));





% segment1 Au
for i=1:Lx1
GnnL=inv((E+1i*eta)*(1+Z0)*eye(4*Lyw)-H_onBDG(tw,uw_Au,Delta0,Lyw,Alphay,(1+Z0)*V_Au)-Vc_BDG'*GnnL*Vc_BDG);
end

% transition region
for i=1:Lx2
H_u=kron(tauz,-1*kron(diag([u(i)*ones(1,Lyw1+Lx2-i),fliplr(u(1:i)),uw_Au*ones(1,Lyw2)]),eye(2)));
H_z=kron(tau0,kron(diag([V_u(i)*ones(1,Lyw1+Lx2-i),fliplr(V_u(1:i)),V_Au*ones(1,Lyw2)]),sgz));
    uw=0;
GnnL=inv((E+1i*eta)*(1+Z0)*eye(4*Lyw)-H_onBDG(tw,uw,Delta0,Lyw,Alphay,0)-H_u-(1+Z0)*H_z-Vc_BDG'*GnnL*Vc_BDG);
end

% segment2 EuS
for i=1:Lx3
H_u=kron(tauz,-1*kron(diag([uw_EuS*ones(1,Lyw1),fliplr(u(1:Lx2)),...
        uw_Au*ones(1,Lyw2)]),eye(2)));
H_z=kron(tau0,kron(diag([V_x1*ones(1,L_y1),V_EuS*ones(1,Lyw1-L_y1),fliplr(V_u(1:Lx2)),...
        V_Au*ones(1,Lyw2)]),sgz));
uw=0;
GnnL=inv((E+1i*eta)*(1+Z0)*eye(4*Lyw)-H_onBDG(tw,uw,Delta0,Lyw,Alphay,0)-H_u-(1+Z0)*H_z-Vc_BDG'*GnnL*Vc_BDG);
end


% segment3 Au
for i=1:Lx1
GnnR=inv((E+1i*eta)*(1+Z0)*eye(4*Lyw)-H_onBDG(tw,uw_Au,Delta0,Lyw,Alphay,(1+Z0)*V_Au)-Vc_BDG*GnnR*Vc_BDG');
end

%transition region



for i=1:Lx2
H_u=kron(tauz,-1*kron(diag([u(i)*ones(1,Lyw1+Lx2-i),fliplr(u(1:i)),uw_Au*ones(1,Lyw2)]),eye(2)));
H_z=kron(tau0,kron(diag([V_u(i)*ones(1,Lyw1+Lx2-i),fliplr(V_u(1:i)),V_Au*ones(1,Lyw2)]),sgz));
uw=0;
GnnR=inv((E+1i*eta)*(1+Z0)*eye(4*Lyw)-H_onBDG(tw,uw,Delta0,Lyw,Alphay,0)-H_u-(1+Z0)*H_z-Vc_BDG*GnnR*Vc_BDG');
 
end

H_u=-1*kron(diag([uw_EuS*ones(1,Lyw1),fliplr(u(1:Lx2)),...
        uw_Au*ones(1,Lyw2)]),eye(2));
H_z=kron(tau0,kron(diag([V_x1*ones(1,L_y1),V_EuS*ones(1,Lyw1-L_y1),fliplr(V_u(1:Lx2)),...
        V_Au*ones(1,Lyw2)]),sgz));

H_onBDG_end=(1+Z0)*H_z+kron(tauz,kron(kron(4*tw*eye(L_leadx),eye(Lyw)),eye(2))+H_u+kron(eye(L_leadx),...
    kron(diag(ones(Lyw-1,1),1),-tw*eye(2)+Alphay*(-1i/2)*sgz)+...
    kron(diag(ones(Lyw-1,1),-1),-tw*eye(2)+Alphay*(1i/2)*sgz))+...
    kron(diag(ones(L_leadx-1,1),1),kron(eye(Lyw),-tw*eye(2)-Alphax*(-1i/2)*sgy))+...
    kron(diag(ones(L_leadx-1,1),-1),kron(eye(Lyw),-tw*eye(2)-Alphax*(1i/2)*sgy)))+...
    kron(taux,kron(eye(Lyw*L_leadx),(Delta0)*eye(2)));
McR=zeros(L_leadx,1);
McR(L_leadx,1)=1;
McL=zeros(1,L_leadx);
McL(1,1)=1;
Vc_BDG_end_R=kron(tauz,kron(McR,kron(eye(Lyw),-tw*eye(2)-Alphax*(-1i/2)*sgy)));
Vc_BDG_end_L=kron(tauz,kron(McL,kron(eye(Lyw),-tw*eye(2)-Alphax*(-1i/2)*sgy)));
gw=inv((E+1i*eta)*(1+Z0)*eye(4*L_leadx*Lyw)-H_onBDG_end-Vc_BDG_end_R*GnnR*Vc_BDG_end_R'-Vc_BDG_end_L'*GnnL*Vc_BDG_end_L);

end