
function Gc_T=tunnel_Au_EuS(k_bT,Lx1,Lx2,Lx3, Lyw1,Lyw2,Lyw,L_y1,uw_Au,uw_EuS,V_x1,V_EuS,V_Au,tw,Alphax,Alphay,Delta,E,spa,D,DD,L_leadx,L_leady,D1,Gamma)



%------------------------------------------------------------------
%position of stm lead
 tauz=[1,0;0,-1];
%-------------------------------------------------------------------
%STM lead
%-------------------------------------------------------------------
tL=3*tw; uL=2*tL;
  tc=tL/5;tc_G=tL/3;
 eta=10^(-7);
HL_BDG=kron(tauz,kron(eye(L_leadx*L_leady),(2*tL-uL)*eye(2))+kron(eye(L_leadx),kron(diag(ones(L_leady-1,1),1)+diag(ones(L_leady-1,1),-1),-tL*eye(2)))+...
   +kron(diag(ones(L_leadx-1,1),1)+diag(ones(L_leadx-1,1),-1),kron(eye(L_leady),-tL*eye(2))));
VcL_BDG=kron(tauz,kron(-tL*eye(L_leadx*L_leady),eye(2)));

Couple1=zeros(L_leady,Lyw);

Couple1(1:L_leady,5+1:5+L_leady)=eye(L_leady); %couple first L_leady sites in the y direction
Couple=kron(eye(L_leadx),Couple1);
Vc_L_S=kron(tauz,kron(-tc*Couple,eye(2)));
%-------------------------------------------------------------------


%-------------------------------------
X=[-2*D:1/1000*DD:-D1,-D1:spa:D1,D1:1/1000*DD:2*D]; 
%X=[0:1/100*DD:0.01];
Sp=[1/1000*DD*ones(1,length(-2*D:1/1000*DD:-D1)),spa*ones(1,length(-D1:spa:D1)),1/1000*DD*ones(1,length(D1:1/1000*DD:2*D))];
%parpool('local',28)
parfor j=1:length(X)
    j
    tic
        %% the grounded lead
%----------------------------------------------------------
HL_BDG_G=kron(tauz,kron(eye(Lyw),(4*tL-2*tL)*eye(2)));
VcL_BDG_G=kron(tauz,kron(eye(Lyw),(-tL)*eye(2)));
G_0n=surleadsts(X(j),eta,HL_BDG_G,VcL_BDG_G);
%----------------------------------------------------------
gL=surleadsts(X(j),eta,HL_BDG,VcL_BDG);
gw=Andreev_end(X(j),eta,tw,uw_Au,uw_EuS,V_x1,V_EuS,V_Au,Alphax,Alphay,Delta,Lx1,Lx2,Lx3, Lyw1,Lyw2,Lyw,L_y1,G_0n,tc_G,L_leadx,Gamma);
Sgma_R=Vc_L_S*gw*Vc_L_S';
Sgma_L=VcL_BDG'*gL*VcL_BDG;
G_F=inv((X(j)+eta*1i)*eye(4*L_leadx*L_leady)-Sgma_R-Sgma_L);
Tau=1i*(Sgma_L-Sgma_L');
r=-eye(4*L_leadx*L_leady)+1i*sqrtm(Tau)*G_F*sqrtm(Tau);
r_ee=r(1:2*L_leadx*L_leady,1:2*L_leadx*L_leady);
r_he=r(2*L_leadx*L_leady+1:4*L_leadx*L_leady,1:2*L_leadx*L_leady);

Gc(j)=abs(trace(eye(2*L_leadx*L_leady)-r_ee'*r_ee+r_he'*r_he));
% 
% G_he(j)=trace(r_he'*r_he);
% G_ee(j)=trace(r_ee'*r_ee);
    toc
end
figure;
plot(X,Gc);
%%--------------------------------------------------------

for j=1:length(E)
W_f=(1/k_bT*exp(-abs(X-E(j))/k_bT))./((exp(-abs((X-E(j)))/k_bT)+1).^2);
A=sum(W_f.*Sp)
Gc_T(1,j)=sum(Gc.*W_f.*Sp)/A;
end
save('long_square_Zeeman_2_point_2Delta.mat');

end


 