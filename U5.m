clear all;
close all;

% Fall 1
% t = 0.0;          % tid (s)
% b = 1.0;          % hjulbasen (m)
% aR = 3.0;         % acceleration höger hjul (m/s^2)
% aL = 3.0;         % acceleration vänster hjul (m/s^2)
% wL = 1.0;         % hastighet vänster hjul (m/s)
% wR = 5.0;         % hastighet höger hjul (m/s)
% x0 = 0.0;         % start x-kordinat
% y0 = 0.0;         % Start y-kordinat
% theta0 = 0.0;     % startvinkel (rad)

% %Fall 2
% t = 0.0;          % tid (s)
% b = 1.0;          % hjulbasen (m)
% aR = 4.0;         % acceleration höger hjul (m/s^2)
% aL = 2.0;         % acceleration vänster hjul (m/s^2)
% wL = 2.0;         % hastighet vänster hjul (m/s)
% wR = 3.0;         % hastighet höger hjul (m/s)
% x0 = 0.0;         % start x-kordinat
% y0 = 0.0;         % Start y-kordinat
% theta0 = 0.0;     % startvinkel (rad)

%Fall 3
% t = 0.0;          % tid (s)
% b = 1.0;          % hjulbasen (m)
% aR = 1.0;         % acceleration höger hjul (m/s^2)
% aL = 3.0;         % acceleration vänster hjul (m/s^2)
% wL = 1.0;         % hastighet vänster hjul (m/s)
% wR = 6.0;         % hastighet höger hjul (m/s)
% x0 = 0.0;         % start x-kordinat
% y0 = 0.0;         % Start y-kordinat
% theta0 = 0.0;     % startvinkel (rad)


%Fall 4 (för uppg. 6)
t = 0.0;          % tid (s)
b = 1.0;          % hjulbasen (m)
aR = 1.2;         % acceleration höger hjul (m/s^2)
aL = 0.4;         % acceleration vänster hjul (m/s^2)
wL = 1.0;         % hastighet vänster hjul (m/s)
wR = 3.0;         % hastighet höger hjul (m/s)
x0 = 0.0;         % start x-kordinat
y0 = 0.0;         % Start y-kordinat
theta0 = 0.0;     % startvinkel (rad)

tslut = 10; % 10 sekunder

function ds = fvel(t, s, b, aL, aR, wL, wR, theta0)
    B = (wR + wL) ./ 2;
    D = (wR - wL) ./ b;
    C = (aR - aL) ./ (2*b);
    A = (aR + aL) ./ 2;

    dxdt = (A * t + B) * cos(C * (t^2) + D * t + theta0);
    dydt = (A * t + B) * sin(C * (t^2) + D * t + theta0);
    dthetadt = 2 * C * t + D;

    ds = [dxdt,dydt,dthetadt]';
end

% test

tspann = [0, tslut];
s0 = [x0, y0, theta0];
options = odeset( 'RelTol' ,1e-6 , 'AbsTol' ,1e-6);
[t,s] = ode45(@(t,s) fvel(t,s,b,aL,aR,wL,wR,theta0), tspann, s0, options);

plot(s(:,1), s(:,2), 'r');
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

%plot(t, [felvektory,felvektorx]);
%grid on
%xlabel('t')
%ylabel('x,y')

%length(t)