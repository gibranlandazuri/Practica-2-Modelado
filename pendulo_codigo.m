function sol = resolverSistema2o()

%parametros
Ip = 0.0079;
Mc = 0.7031;
lp = 0.3302;
Mp = 0.23;
Fc = 0;
Beq = 4.3;
g = 9.81;
Bp = 0.0024;


f = @(t,Y) deal( [] ); %mantener estructura

f = @(t,Y) [ ...
    Y(2); ...
    ( (Ip + Mp*lp^2)*Fc + Mp^2*lp^2*g*cos(Y(3))*sin(Y(3)) - (Ip + Mp*lp^2)*Beq*Y(2) + (Ip*Mp*lp - Mp^2*lp^3)*Y(4)^2*sin(Y(3)) - Mp*lp*Y(4)*cos(Y(3))*Bp ) ...
        / ( (Mc + Mp)*Ip + Mc*Mp*lp^2 + Mp^2*lp^2*sin(Y(3))^2 ); ...
    Y(4); ...
    ( (Mc + Mp) * Mp * g * lp * sin(Y(3)) - (Mc + Mp) * Bp * Y(4) + Fc * Mp * lp * cos(Y(3)) - Mp^2 *lp^2 * Y(4)^2 * sin(Y(3)) * cos(Y(3)) - Beq * Mp * lp * Y(2) * cos(Y(3)) ) ...
        / ( (Mc + Mp) * Ip + Mc * Mp * lp^2 + Mp^2 * lp^2 * sin(Y(3))^2 ) ...
    ];

%condiciones
y1_0 = 0; y1dot_0 = 0; %x dx
y2_0 = 1; y2dot_0 = 0; %alpha dalpha
Y0 = [y1_0; y1dot_0; y2_0; y2dot_0];

tspan = [0 10];


opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

%resolver
[t, Y] = ode45(f, tspan, Y0, opts);

%salida
sol.t = t;
sol.Y = Y; %columnas: [y1, y1', y2, y2']

%graficar
figure;
subplot(2,1,1);
plot(t, Y(:,1), '-b');
xlabel('t'); ylabel('y'); legend('x'); title('Desplazamientos');

subplot(2,1,2);
plot( t, Y(:,3), '-r');
xlabel('t'); ylabel('y'); legend('alpha'); title('Desplazamientos');

figure;
subplot(2,1,1);
plot(t, Y(:,2), '-b');
xlabel('t'); ylabel('y'''); legend('derivada de x'); title('Velocidades');

subplot(2,1,2);
plot( t, Y(:,4), '-r');
xlabel('t'); ylabel('y'''); legend('derivada de alpha'); title('Velocidades');

end