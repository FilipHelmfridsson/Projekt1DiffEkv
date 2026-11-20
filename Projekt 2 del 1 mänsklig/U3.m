function [x, w] = varme_euler(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn)

    % 1. Diskretisering av rummet och tiden
    dx = 1/N; % steglängd i x
    dt = T/M; % tidssteg
    x = (0:N)'; % index
    x = x*dx; % nätpunkter
    j = 2:N; % inre punkter

    w = uinitfcn(x);

    for n = 0:M-1
        t = n*dt; %nuvarande tid

        % 2a, högerledet F(t_n, w_n). så rumsdiskretisering av D*u_xx + f
        uxx = (w(j+1) - 2*w(j) + w(j-1)) / dx^2;   %central differens för u_xx
        F = D*uxx + ffcn(x(j), t);

        % 2b. Euler framåt:
        w(j) = w(j) + dt*F;

        % 3. lite randvillkor
        w(1) = TLfcn(t + dt); % vänster rand, x = 0
        w(end) = TRfcn(t + dt); % höger rand: x = 1
    end
end

%c)
uanalytiskfcn = @(x,t) exp(-25*pi^2*D*t) * sin(5*pi*x);

N = 100; T = 0.01; 
M = T/1e-5;
D = 1;
ffcn = @(x,t) 0*x;
uinitfcn = @(x) sin(5*pi*x);
TLfcn = @(t) 0;
TRfcn = @(t) 0;

[x, w] = varme_euler(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn);

plot(x, abs(w - uanalytiskfcn(x, T)))
xlabel('x'), ylabel('w - u'), title('Fel vid T = 0.01')

maxfel = max(abs(w - uanalytiskfcn(x, T)))


% d)
D = 1; 
M = 1000; % många tidssteg
for N = [20 50 100]
    dx = 1/N;
    dtmax = dx^2/(2*D); % teoretisk gräns
    for faktor = [0.99 1.01 1.05]
        dt = faktor*dtmax;
        T  = M*dt; % sluttid som ger exakt M steg
        [x, w] = varme_euler(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn);
        fprintf('N=%3d, dt=%.2f*dtmax: max|w| = %.2e\n', N, faktor, max(abs(w)))
    end
end