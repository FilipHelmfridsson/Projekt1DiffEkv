clear all;
close all;
clc;

t = 2.0;          % tid (s)
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

%for i = 1:max_iteration
while 1
    f_varde = f(aR)
    df_varde = df(aR)

    aR_ny = aR - f_varde / df_varde;
    fel = abs(aR_ny - aR)
    
    disp(i)
    disp(aR_ny)
    disp(fel)

    if fel < tolerans
        disp('hittat aR!')
        disp( aR_ny)
        break
    end

    aR = aR_ny;
end