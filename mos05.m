clear all; clc; close all;
%% Dane
L1=2;
m1=1; 
m2=2; 
C1=0.1;
R1=10;
k2=1;
k3=2;
b1=4;
b2=1;
b3=2;
l0=1;
d=1;
%% Czas symulacji
t_start = 0;
t_koniec =20;
dt = 0.01;
t = [t_start : dt : t_koniec ];
%% Warunki poczatkowe
x0 = [0 0 0 0 0 0 0 0]';
%% Rozwiazanie dla systemu nieliniowego
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
[t_out, x] = ode23s(@mos5_Kuba_nl, t, x0, options, [L1 m1 m2 C1 R1 k2 k3 b1 b2 b3 l0 d]);

figure(4);
plot(t_out, x(:, :));
legend('£adunek q1 [C]','Natê¿enie pr¹du q1[A]','£adunek q2 [C]','Natê¿enie pr¹du q2[A]','Przemieszczenie m1 [m]', 'Prêdkoœæ m1 [m/s]','Przemieszczenie m2 [m]', 'Prêdkoœæ m2 [m/s]');
title('OdpowiedŸ uk³adu nieliniowego na wymuszenie skokowe')
xlabel('Czas [s]')

figure(6);
subplot(2, 1, 1);
plot(t_out, x(:, :));
legend('£adunek q1 [C]','Natê¿enie pr¹du q1[A]','£adunek q2 [C]','Natê¿enie pr¹du q2[A]','Przemieszczenie m1 [m]', 'Prêdkoœæ m1 [m/s]','Przemieszczenie m2 [m]', 'Prêdkoœæ m2 [m/s]');
title('OdpowiedŸ uk³adu nieliniowego na wymuszenie skokowe')

%% Rozwiazanie dla systemu liniowego
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
[t_out, x] = ode23s(@mos5_Kuba_l, t, x0, options, [L1 m1 m2 C1 R1 k2 k3 b1 b2 b3 l0 d]);

figure(5);
plot(t_out, x(:, :));
legend('£adunek q1 [C]','Natê¿enie pr¹du q1[A]','£adunek q2 [C]','Natê¿enie pr¹du q2[A]','Przemieszczenie m1 [m]', 'Prêdkoœæ m1 [m/s]','Przemieszczenie m2 [m]', 'Prêdkoœæ m2 [m/s]');
title('OdpowiedŸ uk³adu liniowego na wymuszenie skokowe')
xlabel('Czas [s]')
legend('£adunek q1 [C]','Natê¿enie pr¹du q1[A]','£adunek q2 [C]','Natê¿enie pr¹du q2[A]','Przemieszczenie m1 [m]', 'Prêdkoœæ m1 [m/s]','Przemieszczenie m2 [m]', 'Prêdkoœæ m2 [m/s]');
figure(6);
subplot(2, 1, 2);
plot(t_out, x(:, :));
legend('£adunek q1 [C]','Natê¿enie pr¹du q1[A]','£adunek q2 [C]','Natê¿enie pr¹du q2[A]','Przemieszczenie m1 [m]', 'Prêdkoœæ m1 [m/s]','Przemieszczenie m2 [m]', 'Prêdkoœæ m2 [m/s]');
title('OdpowiedŸ uk³adu liniowego na wymuszenie skokowe')
