clear all;clc;
a=10; % lattice constant a amstrong
Delta=0.5;%pairing potential
Gamma=3*Delta;
Z=1/(1+Gamma/Delta); %Renormalization factor from the self energy, at the low energy expansion
tw=Z*1.6*10^4/a^2; %hoping of gold wire
Alphax=Z*0.4*10^3/a; 
Alphay=Z*0.4*10^3/a;
Delta_r=(1-Z)*Delta;
Vz=Delta; %Zeeman gap
tauz=sparse([1,0;0,-1]);
taux=sparse([0,1;1,0]);
tauy=sparse([0,-1i;1i,0]);
sgy=sparse([0,-1i;1i,0]);
sgx=sparse([0,1;1,0]);
sgz=sparse([1,0;0,-1]);
tau0=speye(2);

Lyw=1200/(a); % width of gold wire
Lx1=5000/(a); % lengthth of left side of EuS covered region
Lx3=600/(a); % length of square EuS covered region in oblong  geometry
%Lx2=Lyw-2*Lyw1; % length of intermidiate region between the EuS and gold surface

Ly_EuS=600/a;% width of EuS island


Lx2=0*200/a; % length of transition region
Lyw2=600/a;  %width of bare gold region upper half
Lyw1=Lyw-Ly_EuS-Lyw2;        % width of half segment bare gold surface
  


Lxw=2*Lx2+Lx1*2+Lx3; % total length

V_Au=0.1*Vz; % Zeeman energy of Au region
V_EuS=2*Vz; %Zeeman energy of V_EuS region
uw_Au=200; % Fix the Fermi enery of bare gold region

uw_EuS=26; % tuning the chemical potential of EuS coverd region
hon=zeros(Lxw,Lxw); 
%%
%parpool('local',28)

for j=1:length(uw_EuS)
    j
    tic
    H_N=Hon_squarem(uw_Au,uw_EuS(j),Lx1,Lx2,Lx3,Lyw1,Lyw2,Ly_EuS); % chemecal potential configuration
    H_u=kron(tauz,H_N);
    H_z=Hon_Vz_squarem(V_Au,V_EuS,Lx1,Lx2,Lx3,Lyw1,Lyw2,Ly_EuS); % Zeeman energy configuration
    HBDG=kron(tauz,kron(speye(Lxw),(4*tw)*kron(speye(Lyw),eye(2))+...
        kron(diag(ones(Lyw-1,1),1),-tw*speye(2)+Alphay*(-1i/2)*sgz)+...
        kron(diag(ones(Lyw-1,1),-1),-tw*speye(2)+Alphay*(1i/2)*sgz))+...
        kron(diag(ones(Lxw-1,1),1),kron(speye(Lyw),-tw*eye(2)-Alphax*(-1i/2)*sgy))+...
        kron(diag(ones(Lxw-1,1),-1),kron(speye(Lyw),-tw*eye(2)-Alphax*(1i/2)*sgy)))+...
        Z*H_u+kron(tau0,H_z)+kron(taux,kron(speye(Lyw*Lxw),(Delta_r)*eye(2)));
[V,D]=eigs(HBDG,10,0);
  %E(:,j)=eigs(HBDG,4,0);

    toc
end
%save('E.mat');
%%
% figure;
% plot(uw_EuS,sort(E)/Delta)
% set(gca,'fontsize',20)
% xlabel('\mu(meV)');
% ylabel('E/\Delta');


figure;
set(gca,'fontsize',20)
hold on
Psi=V(:,1).*conj(V(:,1));

psi=zeros(Lyw,Lxw);
for i=1:Lxw
    for j=1:Lyw
        psi(j,i)=Psi((i-1)*Lyw*2+(j-1)*2+1,1)+Psi((i-1)*Lyw*2+(j-1)*2+2,1)+...
            Psi(2*Lyw*Lxw+(i-1)*Lyw*2+(j-1)*2+1,1)+Psi(2*Lyw*Lxw+(i-1)*Lyw*2+(j-1)*2+2,1);
    end
end
[X,Y]=meshgrid(1:Lxw,1:Lyw);
hold on
pcolor(X,Y,psi)
colormap jet
% 
shading interp
% hold on;plot([Lx1+Lx2,Lx1+Lx2],[Lyw1,Lyw1+Ly_EuS],'-r')
% hold on;plot([Lx1+Lx2,Lx1+Lx2+Lx3],[Lyw1+Ly_EuS,Lyw1+Ly_EuS],'-r')
% hold on;plot([Lx1+Lx2+Lx3,Lx1+Lx2+Lx3],[Lyw1+Ly_EuS,Lyw1],'-r')
% hold on;plot([Lx1+Lx2+Lx3,Lx1+Lx2],[Lyw1,Lyw1],'-r')

set(gca,'fontsize',20)
xlabel('Lx');
ylabel('Ly');
title(['Vx=',num2str(Vz),'meV','\Delta=',num2str(Delta),'meV']);
