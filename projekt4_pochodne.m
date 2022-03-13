clear all; clc; close all;
syms z1 dz1 ddz1 z2 dz2 ddz2 q1 dq1 ddq1 q2 dq2 ddq2 L1 L2 m1 m2 C1 R1 k2 k3 b1 b2 b3 u l0 d x1 x2 x3 x4 x5 x6 x7 x8 up x1zero x2zero q1zero q2zero;
K=0.5*L1*dq1^2+0.5*m1*dz1^2+0.5*m2*dz2^2+0.5*((l0/(z1+d))*dq2^2);
U=(1/(2*C1))*q2^2+0.5*k2*(z2-z1)^2+0.5*k3*z2^2;
D=0.5*R1*(dq1-dq2)^2+0.5*b2*(dz2-dz1)^2+0.5*b3*dz2^2+0.5*b1*dz1^2;
L=K-U;

Ldz1=diff(L,dz1);
Ldz1_dt=m1*ddz1;
Lz1=diff(L,z1);
Ddz1=diff(D,dz1);
row1=Ldz1_dt-Lz1+Ddz1;

Ldz2=diff(L,dz2);
Ldz2_dt=m2*ddz2;
Lz2=diff(L,z2);
Ddz2=diff(D,dz2);
row2=Ldz2_dt-Lz2+Ddz2;

Ldq1=diff(L,dq1);
Ldq1_dt=L1*ddq1;
Lq1=diff(L,q1);
Ddq1=diff(D,dq1);
row3=Ldq1_dt-Lq1+Ddq1;

Ldq2=diff(L,dq2)
Ldq2_dt=(l0*(ddq2*(z1+d)-dz1*dq2))/(z1+d)^2;
Lq2=diff(L,q2);
Ddq2=diff(D,dq2);
row4=Ldq2_dt-Lq2+Ddq2;

[sol_ddz1 sol_ddz2 sol_ddq1 sol_ddq2]=solve([row1==0,row2==0,row3==u,row4==0],[ddz1,ddz2,ddq1,ddq2])

%% liniowe
syms z1 dz1 ddz1 z2 dz2 ddz2 q1 dq1 ddq1 q2 dq2 ddq2 L1 L2 m1 m2 C1 R1 k2 k3 b1 b2 b3 u l0 d x1 x2 x3 x4 x5 x6 x7 x8 up x1zero x2zero q1zero q2zero;
x=[x1,x2,x3,x4,x5,x6,x7,x8];
u=[u];

f(1) =x(2);
f(2) =-(x(8)^2*l0 + 2*k2*x(1)^3 + 4*d*k2*x(1)^2 + 2*d^2*k2*x(1) - 2*d^2*k2*x(3) - 2*k2*x(1)^2*x(3) + 2*b1*d^2*x(2) + 2*b2*d^2*x(2) - 2*b2*d^2*x(4) + 2*b1*x(2)*x(1)^2 + 2*b2*x(2)*x(1)^2 - 2*b2*x(4)*x(1)^2 + 4*b1*d*x(2)*x(1) + 4*b2*d*x(2)*x(1) - 4*b2*d*x(4)*x(1) - 4*d*k2*x(1)*x(3))/(2*m1*(d + x(1))^2)
f(3)=x(4);
f(4)=-(b2*x(4) - b2*x(2) + b3*x(4) - k2*x(1) + k2*x(3) + k3*x(3))/m2;
f(5) =x(6);
f(6) =(u - R1*x(6) + R1*x(8))/L1;
f(7)=x(8);
f(8)=-(d^2*x(7) + x(7)*x(1)^2 + 2*d*x(7)*x(1) - C1*R1*d^2*x(6) + C1*R1*d^2*x(8) - C1*R1*x(6)*x(1)^2 + C1*R1*x(8)*x(1)^2 - C1*x(8)*x(2)*l0 - 2*C1*R1*d*x(6)*x(1) + 2*C1*R1*d*x(8)*x(1))/(C1*l0*(d + x(1)));


for i=1:8
    for j=1:8
        A(i,j)=diff(f(i),x(j));
    end
end

for i=1:8
    for j=1:1
          B(i,j)=diff(f(i),u(1));
    end
end

% x1zero=0;
% x2zero=0;
% q1zero=12;
% q2zero=1;
xp=[x1zero 0 x2zero 0 q1zero 0 q2zero 0];
up=0;

A=subs(A,x,xp)
A=subs(A,u,up)
B=subs(B,x,xp)
B=subs(B,u,up)
X=A*(x.'-xp.')+B*(u.'-up.')


