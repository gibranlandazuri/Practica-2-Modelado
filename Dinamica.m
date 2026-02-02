function dx = Dinamica(t, x)

Ip = 0.0079; %parametros
Mc = 0.7031;
lp = 0.3302;
Mp = 0.23;
Fc = 0;
Beq = 4.3;
g = 9.81;
Bp = 0.0024;

xc = x(1); %estados
xc_dot = x(2);
alpha = x(3);
alpha_dot = x(4);

%denominador
D = (Mc+Mp)*Ip + Mc*Mp*lp^2 + Mp^2*lp^2*sin(alpha)^2;

%ecuaciones
xc_ddot = ((Ip + Mp*lp^2)*Fc + Mp^2*lp^2*g*cos(alpha)*sin(alpha) ...
    - (Ip + Mp*lp^2)*Beq*xc_dot - (Ip*Mp*lp - Mp^2*lp^3)*alpha_dot^2*sin(alpha) ...
    - Mp*lp*cos(alpha)*Bp*alpha_dot) / D;

alpha_ddot = ((Mc+Mp)*Mp*g*lp*sin(alpha) - (Mc+Mp)*Bp*alpha_dot + Fc*Mp*lp*cos(alpha) ...
    - Mp^2*lp^2*alpha_dot^2*sin(alpha)*cos(alpha) - Beq*Mp*lp*xc_dot*cos(alpha)) / D;

dx = zeros(4,1); %sistema
dx(1) = xc_dot;
dx(2) = xc_ddot;
dx(3) = alpha_dot;
dx(4) = alpha_ddot;

end
