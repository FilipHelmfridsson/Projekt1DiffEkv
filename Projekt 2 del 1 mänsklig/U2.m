% U2 - Randvärdesproblemet u'' + 16u = f(x),  0 <= x <= 2*pi
%      med randvillkoren u'(0) = gamma och u(2*pi) = delta
clear all; close all;

d = 5;
gamma = 1;
delta = 1;
f = @(x) d*cos(2*x);

% Analytisk lösning från a)
u = @(x) (delta - d/12)*cos(4*x) + gamma/4*sin(4*x) + d/12*cos(2*x);


%% d) Lös med N = 400 och jämför med analytiska lösningen
N = 400;
[x, w] = losBVP(N, f, gamma, delta);

figure;
plot(x, w, x, u(x), '--')
legend('Numerisk', 'Analytisk')
xlabel('x'); ylabel('u(x)')
title('Lösning, N = 400')

figure;
plot(x, w - u(x))
xlabel('x'); ylabel('w - u')
title('Fel, N = 400')

maxfel = max(abs(w - u(x)))


%% e) Konvergensstudie
Nlista = [50 100 200 400 800];
fel = zeros(1, length(Nlista));

for i = 1:length(Nlista)
    N = Nlista(i);
    [x, w] = losBVP(N, f, gamma, delta);

    r = w(2:N) - u(x(2:N));             % felet i de inre punkterna
    fel(i) = sqrt(sum(r.^2) / (N-1));   % normaliserad 2-norm, formel (11)
end

% Felet minskar med 2^p när h halveras  =>  p = log2(fel_N / fel_2N)
ordning = log2(fel(1:end-1) ./ fel(2:end));

disp('       N         fel')
disp([Nlista' fel'])
disp('Ordning mellan på varandra följande N:')
disp(ordning)


%% f) f(x) = exp(-(x-pi)^2) med två olika randvillkor
g = @(x) exp(-(x-pi).^2);

[x, w1] = losBVP(400, g, 0, 0);   % gamma = 0, delta = 0
[x, w2] = losBVP(400, g, 1, 0);   % gamma = 1, delta = 0

figure;
plot(x, w1, x, w2)
legend('\gamma = 0, \delta = 0', '\gamma = 1, \delta = 0')
xlabel('x'); ylabel('u(x)')
title('f(x) = e^{-(x-\pi)^2}')


%% Lösare (finita differenser)
function [x, w] = losBVP(N, f, gamma, delta)
    h = 2*pi/N;
    x = (0:N)' * h;      % alla punkter x_0 ... x_N
    xi = x(2:N);         % inre punkter x_1 ... x_{N-1}
    n = N - 1;           % antal obekanta

    % (w_{j+1} - 2w_j + w_{j-1})/h^2 + 16 w_j = f_j
    A = zeros(n, n);
    for j = 1:n
        A(j,j) = -2/h^2 + 16;
        if j > 1
            A(j,j-1) = 1/h^2;
        end
        if j < n
            A(j,j+1) = 1/h^2;
        end
    end

    % Första raden: w_0 = (4w_1 - w_2 - 2h*gamma)/3 insatt (från formel (6))
    A(1,1) = -2/(3*h^2) + 16;
    A(1,2) =  2/(3*h^2);

    q = f(xi);
    q(1) = q(1) + 2*gamma/(3*h);    % från w_0
    q(n) = q(n) - delta/h^2;        % w_N = delta flyttas till högerledet

    wi = A \ q;

    w0 = (4*wi(1) - wi(2) - 2*h*gamma) / 3;
    w = [w0; wi; delta];
end
