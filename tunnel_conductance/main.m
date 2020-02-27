
clear all;clc;
a=10; 
L_leadx=1; 
L_leady=1;
Delta=0.5;
Vx=Delta;
tw=1.6*10^4/a^2; %hoping of gold wire 
Alphax=0.4*10^3/a;
Alphay=0.4*10^3/a;
V_Au=0.2*Vx; % Zeeman energy of Au region
V_EuS=2.2*Vx; %Zeeman energy of V_EuS region
V_x1=2.2*Vx; %Zeeman energy of EuS covered boundary

Lyw=1200/(a); % length of gold wire
Lx1=5000/(a); % lengthth of left side of EuS covered region
Lx3=600/(a); % length of EuS covered region
Lx2=50/a; % length of intermidiate region between the EuS and gold surface
Lxw=2*Lx2+Lx1*2+Lx3; % total width
Lyw1=600/a; % width of EuS
Lyw2=Lyw-Lx2-Lyw1; % width of bare gold surface
L_y1=100/a;

k_bT=0.04; %temperature
Gamma=3*Delta; %Coupling strength between the superconductor gold surface

uw_Au=200;% Fermi energy of gold surface state
uw_EuS=25; % Chemical potential of EuS part
D=0.5; % energy cut off 
D1=0.05; % a small cut off that needs higher resolution 
DD=0.5; % tuning the spacing
spa=DD/20001; 
E=[-1.2*D:1/1000*DD:-D1,-D1:spa:D1,D1:1/1000*DD:1.2*D];
parpool('local',28)

X=(-2*DD:spa:2*DD); 
% Gc_Vx=zeros(length(X),length(Vx));
for i=1:1:length(Vx)
    i
tic
Gc_T(:,i)=tunnel_Au_EuS(k_bT,Lx1,Lx2,Lx3, Lyw1,Lyw2,Lyw,L_y1,uw_Au,uw_EuS,V_x1,V_EuS,V_Au,tw,Alphax,Alphay,Delta,E,spa,D,DD,L_leadx,L_leady,D1,Gamma);
toc
%%

end
%%
figure
hold on
plot(E,Gc_T,'b')
save('end_square_Zeeman_2_2Delta.mat');

