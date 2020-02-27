function H_N=Hon_squarem(uw_Au,uw_EuS,Lx1,Lx2,Lx3,Lyw1,Lyw2,Ly_EuS)
Lxw=2*Lx2+Lx1*2+Lx3;
Lyw=Lyw1+Lyw2+Ly_EuS;
Ly_c1=Lyw1-Lx2;
Ly_c2=Lyw2-Lx2;
H_N=sparse(2*Lyw*Lxw,2*Lyw*Lxw);
i=1:Lx2;
u=(1-i/(Lx2))*uw_Au+i/(Lx2)*uw_EuS;
v=fliplr(u);




for i=1:Lx1
    H_N(2*(i-1)*Lyw+1:2*(i-1)*Lyw+2*Lyw,2*(i-1)*Lyw+1:2*(i-1)*Lyw+2*Lyw)=-uw_Au*kron(eye(Lyw),eye(2));
end

for i=1:Lx2
    H_N(2*Lx1*Lyw+2*(i-1)*Lyw+1:2*Lx1*Lyw+2*(i-1)*Lyw+2*Lyw,2*Lx1*Lyw+2*(i-1)*Lyw+1:2*Lx1*Lyw+2*(i-1)*Lyw+2*Lyw)=-1*kron(diag([uw_Au*ones(1,Ly_c1),u(1:i),u(i)*ones(1,Lyw-Ly_c1-Ly_c2-2*i),flip(u(1:i)),uw_Au*ones(1,Ly_c2)]),eye(2));
end

for i=1:Lx3
    H_N(2*(Lx1+Lx2)*Lyw+2*(i-1)*Lyw+1:2*(Lx1+Lx2)*Lyw+2*(i-1)*Lyw+2*Lyw,2*(Lx1+Lx2)*Lyw+2*(i-1)*Lyw+1:2*(Lx1+Lx2)*Lyw+2*(i-1)*Lyw+2*Lyw)=-1*kron(diag([uw_Au*ones(1,Ly_c1),u,uw_EuS*ones(1,Ly_EuS),...
       flip(u), uw_Au*ones(1,Ly_c2)]),eye(2));
end

for i=1:Lx2
    H_N(2*(Lx1+Lx2+Lx3)*Lyw+2*(i-1)*Lyw+1:2*(Lx1+Lx2+Lx3)*Lyw+2*(i-1)*Lyw+2*Lyw,2*(Lx1+Lx2+Lx3)*Lyw+2*(i-1)*Lyw+1:2*(Lx1+Lx2+Lx3)*Lyw+2*(i-1)*Lyw+2*Lyw)=...
    -1*kron(diag([uw_Au*ones(1,Ly_c1),u(1,1:Lx2-i+1),v(i)*ones(1,Ly_EuS+2*i-2),flip(u(1,1:Lx2-i+1)),uw_Au*ones(1,Ly_c2)]),eye(2));
end

for i=1:Lx1
    H_N(2*(Lx1+2*Lx2+Lx3)*Lyw+2*(i-1)*Lyw+1:2*(Lx1+2*Lx2+Lx3)*Lyw+2*(i-1)*Lyw+2*Lyw,2*(Lx1+2*Lx2+Lx3)*Lyw+2*(i-1)*Lyw+1:2*(Lx1+2*Lx2+Lx3)*Lyw+2*(i-1)*Lyw+2*Lyw)=-uw_Au*kron(eye(Lyw),eye(2));
end

end