function [x, w] = varme_crank_nicholson(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn)
    dx = 1/N;
    dt = T/M;
    x = (0:N)' * dx;
    xi = x(2:N);

    e = ones(N-1, 1);
    A = D/dx^2 * spdiags([e -2*e e], -1:1, N-1, N-1);

    I = speye(N-1);    % gles identitetsmatris
    L = I - dt/2*A;    % vänster sida
    R = I + dt/2*A;    % höger sida

    b = @(t) ffcn(xi, t) + D/dx^2 * [TLfcn(t); zeros(N-3,1); TRfcn(t)];

    wi = uinitfcn(xi);           % w^0 i de inre punkterna
    for n = 0:M-1
        t  = n*dt;               % gammal tid t_n
        hl = R*wi + dt/2*(b(t) + b(t+dt));   % högerledet, allt känt
        wi = L \ hl;             % lös L*w^{n+1} = hl
    end

    w = [TLfcn(T); wi; TRfcn(T)];   % hela lösningen, N+1 punkter
end

clear; close all;
D = 1;
uanafcn  = @(x,t) sin(exp(t)*x);           % exakt lösning
uinitfcn = @(x) uanafcn(x, 0);
ffcn     = @(x,t) exp(t)*x.*cos(exp(t)*x) + exp(2*t)*sin(exp(t)*x);
TLfcn    = @(t) 0;
TRfcn    = @(t) sin(exp(t));

N = 100; T = 0.1; dt = 1e-3; M = round(T/dt);
[x, w] = varme_crank_nicholson(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn);

plot(x, w - uanafcn(x, T)), xlabel('x'), ylabel('w - u')
maxfel = max(abs(w - uanafcn(x, T)))


% %% e) Konvergensstudie med dt = dx (M = N), T = 1
% T   = 1;
% Ns  = [10 20 40 80 160 320];
% fel = zeros(size(Ns));

% for k = 1:length(Ns)
%     N = Ns(k);
%     [x, w] = varme_crank_nicholson(N, N, T, D, ffcn, uinitfcn, TLfcn, TRfcn);  % M = N
%     r = w(2:end-1) - uanafcn(x(2:end-1), T);   % fel i inre punkter
%     fel(k) = sqrt(mean(r.^2));                 % norm (11)
% end

% ordning = [NaN, log2(fel(1:end-1) ./ fel(2:end))];
% table(Ns', fel', ordning', 'VariableNames', {'N', 'Fel', 'Ordning'})

% % Loglog-plot med referenslinje h^2
% h = 1 ./ Ns;
% figure
% loglog(h, fel, 'o-', h, fel(1)*(h/h(1)).^2, '--')
% xlabel('\Delta x = \Delta t'), ylabel('||r||_2')
% legend('CN', 'h^2 (referens)', 'Location', 'northwest')

N = 100; T = 1; dt = 1e-3; M = round(T/dt);
[x, w] = varme_crank_nicholson(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn);
plot(x, w)
figure;
N = 100; T = 2; dt = 1e-3; M = round(T/dt);
[x, w] = varme_crank_nicholson(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn);
plot(x, w)
figure;
N = 100; T = 4; dt = 1e-3; M = round(T/dt);
[x, w] = varme_crank_nicholson(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn);
plot(x, w)
figure;
N = 100; T = 8; dt = 1e-3; M = round(T/dt);
[x, w] = varme_crank_nicholson(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn);
plot(x, w)