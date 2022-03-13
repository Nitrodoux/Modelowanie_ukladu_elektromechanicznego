function [ dx ] = mos03_odefun_l( t, x, param );
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
%u=3*sin(0.75*2*pi*t);
%% Oblicz pochodne
dx(1,1) =x(2);
dx(2,1) = (-R1/L1)*x(2) + (R1/L1)*x(4) + u/L1;
dx(3,1) = x(4);
dx(4,1) = ((C1*R1)*x(2)*x(5)^2 + (2*C1*R1*d)*x(2)*x(5) + (C1*R1*d^2)*x(2) - x(1)*x(5)^2 + (-2*d)*x(3)*x(5) + (-d^2)*x(3) + (-C1*R1)*x(4)*x(5)^2 + (-2*C1*R1*d)*x(4)*x(5) + (C1*l0)*x(4)*x(6) + (-C1*R1*d^2)*x(4))/((C1*l0)*x(5) + C1*d*l0);
dx(5,1) = x(6);
dx(6,1) = (-l0)*x(4)^2 + (-2*k2)*x(5)^3 + (- 2*b1 - 2*b2)*x(5)^2*x(6) + (2*k2)*x(5)^2*x(7) + (2*b2)*x(5)^2*x(8) + (-4*d*k2)*x(5)^2 + (- 4*b1*d - 4*b2*d)*x(5)*x(6) + (4*d*k2)*x(5)*x(7) + (4*b2*d)*x(5)*x(8) + (-2*d^2*k2)*x(5) + (- 2*b1*d^2 - 2*b2*d^2)*x(5) + (2*d^2*k2)*x(7) + (2*b2*d^2)*x(7);
dx(7,1) = x(8);
dx(8,1) = (k2/m2)*x(5) + (b2/m2)*x(6) + (-(k2 + k3)/m2)*x(7) + (-(b2 + b3)/m2)*x(8);
end

