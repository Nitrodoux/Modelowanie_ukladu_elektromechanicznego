syms l0 d u t ddq1 ddq2 q1 q2 dq1 dq2 z1 z2 dz1 dz2 ddz1 ddz2 L2 R1  C1 m1 m2 b1 b2 b3 k2 k3;
%L2=l0/(d+z1)
K=(1/2)*(L1*dq1^2+(l0/(z1+d))*dq2^2+m1*dz1^2+m2*dz2^2)
U=(1/2)*(1/C1*q2^2+k2*(z2-z1)^2+k3*z2^2)
D=(1/2)*(R1*(dq1-dq2)^2+b1*dz1^2+b2*(dz2-dz1)^2+b3*dz2^2)
L=K-U
%% pochodne dL/dq_d
L_dq1=diff(L,dq1)
L_dq2=diff(L,dq2)
L_dz1=diff(L,dz1)
L_dz2=diff(L,dz2)
%% pochodne dL/dq
L_q1=diff(L,q1)
L_q2=diff(L,q2)
L_z1=diff(L,z1)
L_z2=diff(L,z2)
%% dD/dq_doth
D_dq1=diff(D,dq1)
D_dq2=diff(D,dq2)
D_dz1=diff(D,dz1)
D_dz2=diff(D,dz2)
%% pochodne d/dth dL/dq_doth
dLdq1=L1*ddq1
dLdq2=(l0*(ddq2*(z1+d)-dz1*dq2))/(z1+d)^2
dLdz1=ddz1*m1
dLdz2=ddz2*m2
%% Pe³ne lagrangany
L1p=dLdq1-L_q1+D_dq1-u
L2p=dLdq2-L_q2+D_dq2
L3p=dLdz1-L_z1+D_dz1
L4p=dLdz2-L_z2+D_dz2
%% Wyznaczenie drugich pochodnych
L11=solve(L1p,ddq1)
L22=solve(L2p,ddq2)
L33=solve(L3p,ddz1)
L44=solve(L4p,ddz2)
%% Uporz¹dkowanie
collect(L11,[q1 dq1 q2 dq2 z1 dz1 z2 dz2])
collect(L22,[q1 dq1 q2 dq2 z1 dz1 z2 dz2])
collect(L33,[q1 dq1 q2 dq2 z1 dz1 z2 dz2])
collect(L44,[q1 dq1 q2 dq2 z1 dz1 z2 dz2])
%% Linearyzacja uk³adu
%okreœlenie wektora zmiennych stanu
syms l0 d u t ddq1 ddq2 q1 q2 dq1 dq2 z1 z2 dz1 dz2 ddz1 ddz2 L1 L2 R1  C1 m1 m2 b1 b2 b3 k2 k3;
syms x1 x2 x3 x4 x5 x6 x7 x8;
x=[x1,x2,x3,x4,x5,x6,x7,x8];
u=[u];
%okreœlenie funkcji zmiennych stanu
f(1) = x(2);
f(2) = (-R1/L1)*x(2) + (R1/L1)*x(4) + u/L1;
f(3) = x(4);
f(4) = ((C1*R1)*x(2)*x(5)^2 + (2*C1*R1*d)*x(2)*x(5) + (C1*R1*d^2)*x(2) - x(1)*x(5)^2 + (-2*d)*x(3)*x(5) + (-d^2)*x(3) + (-C1*R1)*x(4)*x(5)^2 + (-2*C1*R1*d)*x(4)*x(5) + (C1*l0)*x(4)*x(6) + (-C1*R1*d^2)*x(4))/((C1*l0)*x(5) + C1*d*l0);
f(5) = x(6);
f(6) = (-l0)*x(4)^2 + (-2*k2)*x(5)^3 + (- 2*b1 - 2*b2)*x(5)^2*x(6) + (2*k2)*x(5)^2*x(7) + (2*b2)*x(5)^2*x(8) + (-4*d*k2)*x(5)^2 + (- 4*b1*d - 4*b2*d)*x(5)*x(6) + (4*d*k2)*x(5)*x(7) + (4*b2*d)*x(5)*x(8) + (-2*d^2*k2)*x(5) + (- 2*b1*d^2 - 2*b2*d^2)*x(5) + (2*d^2*k2)*x(7) + (2*b2*d^2)*x(7);
f(7) = x(8);
f(8) = (k2/m2)*x(5) + (b2/m2)*x(6) + (-(k2 + k3)/m2)*x(7) + (-(b2 + b3)/m2)*x(8);
%% obliczenie pochodnych cz¹stkowych
for i=1:8
    for j=1:8
        A(i,j)=diff(f(i),x(j));
    end
end
%%
for i=1:8
    for j=1:1
          B(i,j)=diff(f(i),u(1));
    end
end

q1zero=0;
q2zero=0;
x1zero=0;
x2zero=0;
pp=[q1zero 0 q2zero 0 x1zero 0 x2zero 0];
up=0;
%Podstawienie punktów pracy do macierzy A
A=subs(A,x,pp)
A=subs(A,u,up)
B=subs(B,x,pp)
B=subs(B,u,up)
% wyznaczenie rownan liniowych przy pomocy macierzy A
X=A*(x.'-pp.')
X1=B*(u.'-up.')
