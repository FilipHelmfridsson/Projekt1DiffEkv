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

tslut = 7*pi/8;

% Utifrån U2 löste vi ut Tlap som är tiden för ett varv Tslut =Tlap (av våra konstanter i denna uppgift)
% Ta in vektor S, returnera derivata med avseende på tid.

ds = fvel(t ,s ,b , aL , aR , wL , wR, theta0);

%% _____ Euler _____

felvektor = [];

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

for k = [0,1,2]
    heuler = tslut/100;
    heuler = heuler / (2^k);

    s = [x0,y0,theta0]; % stratposition och vinkelB

    for i = [t: heuler: tslut-heuler];
        s(end+1, :) = s(end, :) + heuler * fvel(i, s(end, :), b, aL, aR, wL, wR, theta0);
    end

    %disp(s);
    plot(s(:,1),s(:,2), 'r')
    axis equal;

    xeuler = s(end,1);


    %%________Runge kutta 4________%
    hrk4 = tslut/100;
    hrk4 = hrk4 / (2^k);

    s = [x0,y0,theta0];

    for i = [t: hrk4: tslut-hrk4];
        k1 = fvel(i, s(end, :), b, aL, aR, wL, wR, theta0);

        k2 = fvel(i + hrk4 / 2, s(end, :) + (hrk4 / 2) .* k1, b, aL, aR, wL, wR, theta0);

        k3 = fvel(i + hrk4 / 2, s(end, :) + (hrk4 / 2) .* k2, b, aL, aR, wL, wR, theta0);

        k4 = fvel(i + hrk4, s(end, :) + hrk4 .* k3, b, aL, aR, wL, wR, theta0);
        
        s(end+1, :) = s(end, :) + (hrk4 / 6) .* (k1 +  2 .* k2 + 2 .* k3 + k4);
    end

    %disp(s);
    plot(s(:,1),s(:,2), 'b');
    axis equal;

    xrk4 = s(end,1);


    %radie = abs(b.*(wR+wL)./(2.*(wR-wL)));

    %analytiskx = sin(-tslut) * radie


    % Analytisk metod
    B = (wR + wL) ./ 2;
    D = (wR - wL) ./ b;

    x_exakt = @(t) x0 + (B/D) * (sin(D*t + theta0) - sin(theta0));
    y_exakt = @(t) y0 - (B/D) * (cos(D*t + theta0) - cos(theta0));


    x_slut = x_exakt(tslut);
    y_slut = y_exakt(tslut);

    felrk4 = abs(xrk4 - x_slut)
    feleuler = abs(xeuler - x_slut)

    felvektor = [felvektor;2^k,felrk4,feleuler];
end

noggranhetsordningrk4 = log((abs(felvektor(1,2) - felvektor(2,2))) / abs((felvektor(3,2) - felvektor(2,2))))/log(2)
noggranhetsordningeuler = log((abs(felvektor(1,3) - felvektor(2,3))) / abs((felvektor(3,3) - felvektor(2,3))))/log(2)

p = log(felvektor(1,2) / felvektor(2,2)) / log(2);
p = log(felvektor(1,3) / felvektor(2,3)) / log(2);