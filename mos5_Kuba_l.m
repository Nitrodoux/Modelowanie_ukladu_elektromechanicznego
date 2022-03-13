function [ dx ] = mos03_odefun_l( t, x, param )
%% Pobierz wartosci parametrow
L1=param(1);
m1=param(2); 
m2=param(3); 
C1=param(4);
R1=param(5);
k2=param(6);
k3=param(7);
b1=param(8);
b2=param(9);
b3=param(10);
l0=param(11);
d=param(12);
%% Wyznacz wymuszenie
u=1;
q1zero=0;
q2zero=0;
x1zero=0;
x2zero=0;
%u=3*sin(0.75*2*pi*t);
%% Oblicz pochodne
%{
dx(1,1)=x(2);
dx(2,1)= u/L1 - (R1*x(2))/L1 + (R1*x(4))/L1;
dx(3,1)=x(4);
dx(4,1)=((q2zero - x(3))*(d^2 + 2*d*x1zero + x1zero^2))/(C1*l0*(d + x1zero)) - ((2*d*q2zero + 2*q2zero*x1zero)/(C1*l0*(d + x1zero)) - (q2zero*d^2 + 2*q2zero*d*x1zero + q2zero*x1zero^2)/(C1*l0*(d + x1zero)^2))*(x(5) - x1zero) + (x(2)*(C1*R1*d^2 + 2*C1*R1*d*x1zero + C1*R1*x1zero^2))/(C1*l0*(d + x1zero)) - (x(4)*(C1*R1*d^2 + 2*C1*R1*d*x1zero + C1*R1*x1zero^2))/(C1*l0*(d + x1zero));
dx(5,1)=x(6);
dx(6,1)=(x(5) - x1zero)*((2*k2*d^2*x1zero - 2*k2*x2zero*d^2 + 4*k2*d*x1zero^2 - 4*k2*x2zero*d*x1zero + 2*k2*x1zero^3 - 2*k2*x2zero*x1zero^2)/(m1*(d + x1zero)^3) - (2*k2*d^2 + 8*k2*d*x1zero - 4*k2*x2zero*d + 6*k2*x1zero^2 - 4*k2*x2zero*x1zero)/(2*m1*(d + x1zero)^2)) - (x(6)*(2*b1*d^2 + 2*b2*d^2 + 2*b1*x1zero^2 + 2*b2*x1zero^2 + 4*b1*d*x1zero + 4*b2*d*x1zero))/(2*m1*(d + x1zero)^2) + (x(8)*(2*b2*d^2 + 4*b2*d*x1zero + 2*b2*x1zero^2))/(2*m1*(d + x1zero)^2) + ((x(7) - x2zero)*(2*k2*d^2 + 4*k2*d*x1zero + 2*k2*x1zero^2))/(2*m1*(d + x1zero)^2);
dx(7,1)=x(8);
dx(8,1)= (b2*x(6))/m2 - ((k2 + k3)*(x(7) - x2zero))/m2 + (k2*(x(5) - x1zero))/m2 - (x(8)*(b2 + b3))/m2;
%}

dx(1,1)=x(2);
dx(2,1)= (R1*x(4))/L1 - (R1*x(2))/L1+u/L1;
dx(3,1)=x(4);
dx(4,1)=(R1*d*x(2))/l0 - (R1*d*x(4))/l0 - (d*x(3))/(C1*l0);
dx(5,1)=x(6);
dx(6,1)=x(7)*(2*b2*d^2 + 2*d^2*k2) - x(5)*(2*b1*d^2 + 2*b2*d^2 + 2*d^2*k2);
dx(7,1)=x(8);
dx(8,1)=(b2*x(6))/m2 + (k2*x(5))/m2 - (x(8)*(b2 + b3))/m2 - (x(7)*(k2 + k3))/m2;
%}
end
