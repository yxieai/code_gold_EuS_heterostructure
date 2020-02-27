function  gs=surleadsts(E,eta,Hins,As)

Ds=(E+1i*eta)*eye(length(Hins))-Hins;
Bs=As';
Di=inv(Ds);
d=Ds;
while (sum(sum(As.*conj(As)))>10^(-10))
   d=d-As*Di*Bs;
   Ds=Ds-As*Di*Bs-Bs*Di*As;
   Bs=Bs*Di*Bs; 
   As=As*Di*As;
   Di=inv(Ds);
end
 gs=inv(d);
end