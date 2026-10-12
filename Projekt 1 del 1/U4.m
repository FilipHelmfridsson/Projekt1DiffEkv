clear all;
close all;

t = 0.0;          % tid (s)
b = 1.0;          % hjulbasen (m)
aR = 0.0;         % acceleration höger hjul (m/s^2)
aL = 0.0;         % acceleration vänster hjul (m/s^2)
wL = 4.0;         % hastighet vänster hjul (m/s)
wR = 2.0;         % hastighet höger hjul (m/s)
x0 = 0.0;         % start x-kordinat
y0 = 0.0;         % Start y-kordinat
theta0 = 0.0;     % startvinkel (rad)

tslut = pi;


tspann = [0, tslut];
s0 = [x0, y0, theta0];
options = odeset( 'RelTol', 1e-8, 'AbsTol', 1e-8, 'Refine', 5);
[t,s] = ode45(@(t,s) fvel(t,s,b,aL,aR,wL,wR)', tspann, s0, options);

disp("Antal tidssteg:")
disp(length(t))

plot(s(:,1), s(:,2), 'k')
axis equal

% Analytisk metod
B = (wR + wL) ./ 2;
D = (wR - wL) ./ b;

x_exakt = @(t) x0 + (B/D) * (sin(D*t + theta0) - sin(theta0));
y_exakt = @(t) y0 - (B/D) * (cos(D*t + theta0) - cos(theta0));

xlosningar = x_exakt(t);
ylosningar = y_exakt(t);

felvektorx = abs(xlosningar - s(:,1));
felvektory = abs(ylosningar - s(:,2));

plot(t, [felvektory,felvektorx]);
grid on
xlabel('t')
ylabel('x,y')
