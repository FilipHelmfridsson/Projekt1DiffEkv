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
h = pi/100;

% Utifrån U2 löste vi ut Tlap som är tiden för ett varv Tslut =Tlap (av våra konstanter i denna uppgift)
% Ta in vektor S, returnera derivata med avseende på tid.
s = [x0,y0,theta0]; % stratposition och vinkel

ds = fvel(t ,s ,b , aL , aR , wL , wR, theta0);

%% _____ Euler _____


for i = [t: h: tslut];
    s(end+1, :) = s(end, :) + h * fvel(i, s(end, :), b, aL, aR, wL, wR, theta0);
end

%disp(s);
plot(s(:,1),s(:,2));
axis equal;
    
function ds = fvel(t, s, b, aL, aR, wL, wR, theta0)
    B = (wR + wL) ./ 2;
    D = (wR - wL) ./ b;
    C = (aR - aL) ./ (2*b);
    A = (aR + aL) ./ 2;

    dxdt = (A * t + B) * cos(C * (t^2) + D * t + theta0);
    dydt = (A * t + B) * sin(C * (t^2) + D * t + theta0);
    dthetadt = 2 * C * t + D;

    ds = [dxdt,dydt,dthetadt];
end