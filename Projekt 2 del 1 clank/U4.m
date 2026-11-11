% U4: Värmeledningsekvationen med Crank–Nicolson, konstruerad lösning (17)-(18)
clear; close all;

D = 1;
uanafcn  = @(x,t) sin(exp(t)*x);
uinitfcn = @(x) uanafcn(x, 0);
ffcn     = @(x,t) exp(t)*x.*cos(exp(t)*x) + exp(2*t)*sin(exp(t)*x);
TLfcn    = @(t) 0;
TRfcn    = @(t) sin(exp(t));
solve    = @(N, M, T) heat_cn(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn);

%% c) N = 100, dt = 1e-3, T = 0.1
T = 0.1;
[x, w] = solve(100, round(T/1e-3), T);
plot(x, w - uanafcn(x, T))
xlabel('x'), ylabel('w - u'), title('Fel vid T = 0.1')
maxfel = max(abs(w - uanafcn(x, T)))

%% d) Stabilitet: stora tidssteg, långt över Euler-gränsen dx^2/2 = 5e-5
N = 100; T = 1;
for M = [1 5 10 100]
    [x, w] = solve(N, M, T);
    fprintf('dt = %5.2f:  max fel = %.2e\n', T/M, max(abs(w - uanafcn(x, T))))
end

%% e) Konvergensstudie med dt = dx (M = N), T = 1
Ns  = [10 20 40 80 160 320];
fel = zeros(size(Ns));
for k = 1:length(Ns)
    [x, w] = solve(Ns(k), Ns(k), 1);
    r = w(2:end-1) - uanafcn(x(2:end-1), 1);
    fel(k) = sqrt(mean(r.^2));
end
ordning = [NaN, log2(fel(1:end-1) ./ fel(2:end))];
table(Ns', fel', ordning', 'VariableNames', {'N', 'Fel', 'Ordning'})

%% f) Den konstruerade lösningen vid t = 1, 2, 4
xf = linspace(0, 1, 2000);
figure, hold on
for t = [1 2 4]
    plot(xf, uanafcn(xf, t))
end
legend('t = 1', 't = 2', 't = 4'), xlabel('x'), ylabel('u(x,t)')
