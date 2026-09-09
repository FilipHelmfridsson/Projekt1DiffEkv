b = 1.0;          % Hjulbas (m)
wL = 4.0;         % Hastighet vänster hjul (m/s)
wR = 2.0;         % Hastighet höger hjul (m/s)
x0 = 0.0;         % Start x-koordinat
y0 = 0.0;         % Start y-koordinat
theta0 = 0.0;     % Startvinkel (rad)

[xbar, ybar, R, Tlap] = robotcircle(b, wL, wR, x0, y0, theta0)

function [xbar, ybar, R, Tlap] = robotcircle(b, wL, wR, x0, y0, theta0)

    B = (wR + wL) ./ 2;
    D = (wR - wL) ./ b;


    xbar = x0 - (B / D) * sin(theta0);

    ybar = y0 - (B / D) * cos(theta0);
    
    Tlap = (2. * pi .* b) / abs(wR - wL);
    

    R = abs(B / D);
    
end