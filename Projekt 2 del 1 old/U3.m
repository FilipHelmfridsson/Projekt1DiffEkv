% U3: Värmeledningsekvationen med Euler framåt
clear; close all;

D = 1;
ffcn     = @(x,t) 0*x;
uinitfcn = @(x) sin(5*pi*x);
TLfcn    = @(t) 0;
TRfcn    = @(t) 0;
uanafcn  = @(x,t) exp(-25*pi^2*D*t) * sin(5*pi*x);

%% c) N = 100, dt = 1e-5, T = 0.01
N = 100; T = 0.01; M = round(T/1e-5);
[x, w] = heat_euler(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn);
plot(x, w - uanafcn(x, T))
xlabel('x'), ylabel('w - u'), title('Fel vid T = 0.01')
maxfel = max(abs(w - uanafcn(x, T)))

%% d) Stabilitet: dt_max = dx^2/(2D). Kör M = 1000 steg strax under/över gränsen.
M = 1000;
for N = [20 50 100]
    dtmax = (1/N)^2 / (2*D);
    for faktor = [0.99 1.01 1.05]
        [~, w] = heat_euler(N, M, M*faktor*dtmax, D, ffcn, uinitfcn, TLfcn, TRfcn);
        fprintf('N = %3d, dt = %.2f*dt_max:  max|w| = %.2e\n', N, faktor, max(abs(w)))
    end
end
