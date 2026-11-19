function [x, w] = bvp_fdm(N, f, gamma, delta)
% BVP_FDM  Löser randvärdesproblemet
%
%       u''(x) + 16 u(x) = f(x),     0 <= x <= 2*pi
%       u'(0)  = gamma               (Neumann-villkor i vänster rand)
%       u(2*pi) = delta              (Dirichlet-villkor i höger rand)
%
%   med finita differenser (centrala differenser i de inre punkterna).
%
%   Indata:
%       N      - antal delintervall (ger N+1 gitterpunkter)
%       f      - funktionshandtag för högerledet, t.ex. @(x) sin(x)
%       gamma  - värdet på u'(0)
%       delta  - värdet på u(2*pi)
%
%   Utdata:
%       x  - kolumnvektor med alla N+1 punkter x_0, ..., x_N
%       w  - kolumnvektor med den approximativa lösningen w_j ≈ u(x_j)


%% ------------------------------------------------------------------------
%  Steg 1: Skapa gittret
%  ------------------------------------------------------------------------
%  Intervallet [0, 2*pi] delas in i N lika stora delintervall med steglängd h.

h = 2*pi / N;

x = (0:N)' * h;       % alla punkter:  x_0 = 0, x_1 = h, ..., x_N = 2*pi


%% ------------------------------------------------------------------------
%  Steg 2: Bestäm vilka värden som är okända
%  ------------------------------------------------------------------------
%  - w_N är känd direkt från randvillkoret:  w_N = u(2*pi) = delta.
%  - w_0 är INTE känd, men kan uttryckas med w_1 och w_2 via randvillkoret
%    u'(0) = gamma (se steg 4). Den elimineras alltså ur systemet.
%  Kvar att lösa för blir de inre punkterna w_1, ..., w_{N-1}.

xi = x(2:N);          % inre punkter x_1, ..., x_{N-1}
n  = N - 1;           % antal obekanta = storleken på ekvationssystemet


%% ------------------------------------------------------------------------
%  Steg 3: Bygg matrisen för de inre punkterna
%  ------------------------------------------------------------------------
%  I varje inre punkt x_j ersätts u'' med en central differens:
%
%       u''(x_j) ≈ ( w_{j+1} - 2 w_j + w_{j-1} ) / h^2
%
%  Insatt i differentialekvationen ger det ekvationen
%
%       (1/h^2) w_{j-1} + (-2/h^2 + 16) w_j + (1/h^2) w_{j+1} = f(x_j)
%
%  Varje rad i matrisen A får alltså tre element (tridiagonal matris):
%       underdiagonal:  1/h^2
%       diagonal:      -2/h^2 + 16
%       överdiagonal:   1/h^2

A = zeros(n, n);

for j = 1:n

    % diagonalelementet (koefficient framför w_j)
    A(j, j) = -2/h^2 + 16;

    % elementet till vänster (koefficient framför w_{j-1}),
    % finns inte på första raden
    if j > 1
        A(j, j-1) = 1/h^2;
    end

    % elementet till höger (koefficient framför w_{j+1}),
    % finns inte på sista raden
    if j < n
        A(j, j+1) = 1/h^2;
    end

end


%% ------------------------------------------------------------------------
%  Steg 4: Hantera vänster randvillkor  u'(0) = gamma
%  ------------------------------------------------------------------------
%  Derivatan i x_0 approximeras med en ensidig differens av andra ordningen
%  (formel (6)):
%
%       u'(0) ≈ ( -3 w_0 + 4 w_1 - w_2 ) / (2h) = gamma
%
%  Löser vi ut w_0 får vi
%
%       w_0 = ( 4 w_1 - w_2 - 2 h gamma ) / 3
%
%  Detta sätts in i första ekvationen (j = 1):
%
%       ( w_2 - 2 w_1 + w_0 ) / h^2 + 16 w_1 = f(x_1)
%
%  Efter förenkling blir det
%
%       ( -2/(3h^2) + 16 ) w_1 + ( 2/(3h^2) ) w_2 = f(x_1) + 2 gamma / (3h)
%
%  Första raden i A skrivs därför över med de nya koefficienterna.

A(1, 1) = -2/(3*h^2) + 16;
A(1, 2) =  2/(3*h^2);


%% ------------------------------------------------------------------------
%  Steg 5: Bygg högerledet
%  ------------------------------------------------------------------------
%  Grundidén: högerledet är f utvärderad i de inre punkterna.

q = f(xi);

%  Första raden: termen med gamma från steg 4 hamnar i högerledet.
q(1) = q(1) + 2*gamma / (3*h);

%  Sista raden (j = N-1): ekvationen innehåller w_N, som är känd (= delta).
%  Termen (1/h^2) * delta flyttas över till högerledet.
q(n) = q(n) - delta / h^2;


%% ------------------------------------------------------------------------
%  Steg 6: Lös ekvationssystemet  A * wi = q
%  ------------------------------------------------------------------------

wi = A \ q;           % lösningen i de inre punkterna w_1, ..., w_{N-1}


%% ------------------------------------------------------------------------
%  Steg 7: Sätt ihop hela lösningen
%  ------------------------------------------------------------------------
%  - w_0 räknas ut från formeln i steg 4 (nu när w_1 och w_2 är kända).
%  - w_N är randvärdet delta.

w0 = (4*wi(1) - wi(2) - 2*h*gamma) / 3;

w = [w0; wi; delta];  % alla N+1 värden: w_0, w_1, ..., w_N

end
