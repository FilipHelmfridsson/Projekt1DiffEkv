clear all;
close all;
clc;

t = 1.0;          % tid (s)
b = 1.0;          % hjulbasen (m)
aL = 0.4;         % acceleration vänster hjul (m/s^2)
wL = 1.0;         % hastighet vänster hjul (m/s)
wR = 3.0;         % hastighet höger hjul (m/s)
R = 1.0; % Radien

vL = aL * t + wL;
    
%R = @(t) (b * ( (aR*t + wR) + vL )) ./ (2 * ( (aR*t + wR) - vL ))


f = @(aR) (b * ( (aR*t + wR) + vL )) ./ (2 * ( (aR*t + wR) - vL )) - R;  % Nollställesfunktionen
df = @(aR) - (b * t * vL) ./ ((aR*t + wR - vL).^2 ); % Derivatan för nollställefunktionen


% (1) Plot av aR
aR_test = [0:1:5];

svarsvektor = [f(aR_test)]

plot(aR_test,svarsvektor)

grid on;
hold on;

% (2) Newton

aR = 1.0; % Gissad, vald till 1 utifrån plot.
tolerans = 1e-6;

i = 1;
felGammal = 1;
felGammalGammal = 1;
while 1
    f_varde = f(aR);
    df_varde = df(aR);

    aR_ny = aR - f_varde / df_varde;
    fel = abs(aR_ny - aR);
    
    disp("Iteration:")
    disp(i)
    disp("Nytt aR värde:")
    disp(aR_ny)
    disp("Felet är:")
    disp(fel)
    if i >= 3
        konvergensOrdning = log(fel / felGammal) / log(felGammal / felGammalGammal);
        disp("Konvergensordningen är:")
        disp(konvergensOrdning)
    else
        konvergensOrdning = 0;
        disp("Det krävs ett par iterationer för att beräkna konvergensordning")
    end

    if fel < tolerans
        disp('hittat aR!')
        disp(aR_ny)
        break
    end

    aR = aR_ny;
    felGammalGammal=felGammal;
    felGammal = fel;
    i = i+1;
end