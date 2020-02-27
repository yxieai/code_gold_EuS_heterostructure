function H_NZ=Hon_Vz_squarem(V_Au,V_EuS,Lx1,Lx2,Lx3,Lyw1,Lyw2,Ly_EuS)
Lxw=Lx1*2+Lx3+2*Lx2;
Lyw=Lyw1+Lyw2+Ly_EuS;
H_NZ=sparse(2*Lyw*Lxw,2*Lyw*Lxw);
Ly_c1=Lyw1-Lx2;
Ly_c2=Lyw2-Lx2;
sgz=[1,0;0,-1];
i=1:Lx2;
u=(1-i/(Lx2))*V_Au+i/(Lx2)*V_EuS;
v=fliplr(u);


%segment 1
for i=1:Lx1
    H_NZ(2*(i-1)*Lyw+1:2*(i-1)*Lyw+2*Lyw,2*(i-1)*Lyw+1:2*(i-1)*Lyw+2*Lyw)=V_Au*kron(eye(Lyw),sgz);
end


% transition

for i=1:Lx2
    H_NZ(2*(Lx1)*Lyw+2*(i-1)*Lyw+1:2*(Lx1)*Lyw+2*(i-1)*Lyw+2*Lyw,2*(Lx1)*Lyw+2*(i-1)*Lyw+1:2*(Lx1)*Lyw+2*(i-1)*Lyw+2*Lyw)=kron(diag([V_Au*ones(1,Ly_c1),u(1:i),u(i)*ones(1,Lyw-Ly_c1-Ly_c2-2*i),flip(u(1:i)),V_Au*ones(1,Ly_c2)]),sgz);
end


%segment2 EuS covered gold
for i=1:Lx3
    H_NZ(2*(Lx1+Lx2)*Lyw+2*(i-1)*Lyw+1:2*(Lx1+Lx2)*Lyw+2*(i-1)*Lyw+2*Lyw,2*(Lx1+Lx2)*Lyw+2*(i-1)*Lyw+1:2*(Lx1+Lx2)*Lyw+2*(i-1)*Lyw+2*Lyw)=kron(diag([V_Au*ones(1,Ly_c1),u,V_EuS*ones(1,Ly_EuS),...
       flip(u), V_Au*ones(1,Ly_c2)]),sgz);
end

% transition

for i=1:Lx2
   H_NZ(2*(Lx1+Lx2+Lx3)*Lyw+2*(i-1)*Lyw+1:2*(Lx1+Lx2+Lx3)*Lyw+2*(i-1)*Lyw+2*Lyw,2*(Lx1+Lx2+Lx3)*Lyw+2*(i-1)*Lyw+1:2*(Lx1+Lx2+Lx3)*Lyw+2*(i-1)*Lyw+2*Lyw)=...
   kron(diag([V_Au*ones(1,Ly_c1),u(1,1:Lx2-i+1),v(i)*ones(1,Ly_EuS+2*i-2),flip(u(1,1:Lx2-i+1)),V_Au*ones(1,Ly_c2)]),sgz);
end

%segment 2
for i=1:Lx1
   H_NZ(2*(Lx1+2*Lx2+Lx3)*Lyw+2*(i-1)*Lyw+1:2*(Lx1+2*Lx2+Lx3)*Lyw+2*(i-1)*Lyw+2*Lyw,2*(Lx1+2*Lx2+Lx3)*Lyw+2*(i-1)*Lyw+1:2*(Lx1+2*Lx2+Lx3)*Lyw+2*(i-1)*Lyw+2*Lyw)=V_Au*kron(eye(Lyw),sgz);
end

end