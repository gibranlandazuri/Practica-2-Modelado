clc; clear; close all

T = [0 10];

xc0     = 0;
xc_dot0 = 0;
alpha0  = 1*pi/180;   %radianes
alpha_dot0 = 0;

x0 = [xc0; xc_dot0; alpha0; alpha_dot0];

[t, x] = ode45(@Dinamica, T, x0);

%graficas
figure
subplot(2,1,1)
plot(t, x(:,1), 'LineWidth',1.5)
grid on
xlabel('Tiempo [s]')
ylabel('x_c [m]')
title('Posición del carro')

subplot(2,1,2)
plot(t, x(:,3)*180/pi, 'LineWidth',1.5)
grid on
xlabel('Tiempo [s]')
ylabel('\alpha [deg]')
title('Ángulo del péndulo')
