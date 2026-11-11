% U2: Randvärdesproblemet u'' + 16u = f, u'(0) = gamma, u(2*pi) = delta
clear; close all;

d = 5; gamma = 1; delta = 1;
f    = @(x) d*cos(2*x);
uana = @(x) (delta - d/12)*cos(4*x) + gamma/4*sin(4*x) + d/12*cos(2*x);

%% d) N = 400
[x, w] = bvp_fdm(400, f, gamma, delta);
figure
plot(x, w, x, uana(x), '--')
legend('Numerisk', 'Analytisk'), xlabel('x'), ylabel('u')

figure
plot(x, w - uana(x))
xlabel('x'), ylabel('w - u'), title('Fel, N = 400')

%% e) Konvergensstudie
Nlista = [50 100 200 400 800];
fel = zeros(1, length(Nlista));
for i = 1:length(Nlista)
    N = Nlista(i);
    [x, w] = bvp_fdm(N, f, gamma, delta);
    r = w(2:N) - uana(x(2:N));             % felet i de inre punkterna
    fel(i) = sqrt(sum(r.^2) / (N-1));      % normaliserad 2-norm, formel (11)
end

% Felet minskar med 2^p när h halveras  =>  p = log2(fel_N / fel_2N)
ordning = log2(fel(1:end-1) ./ fel(2:end));

disp('       N         fel')
disp([Nlista' fel'])
disp('Ordning mellan på varandra följande N:')
disp(ordning)

%% f) f = exp(-(x-pi)^2)
g = @(x) exp(-(x - pi).^2);
[x, w1] = bvp_fdm(400, g, 0, 0);
[x, w2] = bvp_fdm(400, g, 1, 0);
figure
plot(x, w1, x, w2)
legend('\gamma = 0, \delta = 0', '\gamma = 1, \delta = 0'), xlabel('x'), ylabel('u')
