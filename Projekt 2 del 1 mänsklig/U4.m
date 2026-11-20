clear all;
close all;

D = 1;
u_analytisk = @(x,t) sin(exp(t)*x); %exakt lösning
begynnelse_varme = @(x) u_analytisk(x, 0);
f = @(x,t) exp(t)*x.*cos(exp(t)*x) + exp(2*t)*sin(exp(t)*x);
t_vanster = @(t) 0;
t_hoger = @(t) sin(exp(t));

N = 100; T = 0.1; dt = 1e-3; M = round(T/dt);
[x, t, U] = varmeledning_theta(D, 1, N, T, M, 0.5, begynnelse_varme, t_vanster, t_hoger, f);
w = U(end,:).'; %Tar sista radien alla klumner för att få stavens temperatur i sista tidssteget

plot(x, w - u_analytisk(x, T)), xlabel('x'), ylabel('w - u')
maxfel = max(abs(w - u_analytisk(x, T)))
figure;


%% e) Konvergensstudie med dt = dx (M = N), T = 1
T   = 1;
Ns  = [10 20 40 80 160 320];
fel = zeros(size(Ns));

for k = 1:length(Ns)
    N = Ns(k);
    [x, t, U] = varmeledning_theta(D, 1, N, T, N, 0.5, begynnelse_varme, t_vanster, t_hoger, f); %M = N
    w = U(end,:).';
    r = w(2:end-1) - u_analytisk(x(2:end-1), T); %fel i inre punkter
    fel(k) = sqrt(mean(r.^2)); %norm (11)
end

ordning = [0, log2(fel(1:end-1) ./ fel(2:end))]; %Lägg till nolla i början då det inte går at tberäkna ordningen för första mätningen
for i = 1:length(Ns)
    fprintf('För Ns = %d är felet %.3e och ordningen %.2f.\n', Ns(i), fel(i), ordning(i));
end

%Loglog-plot med referenslinje h^2
h = 1 ./ Ns;

loglog(h, fel, h, fel(1)*(h/h(1)).^2)
xlabel('Delta x = Delta t'), ylabel('||r||_2')
legend('CN', 'h^2 (referens)')
axis equal;
figure;

N = 100; T = 1; dt = 1e-3; M = round(T/dt);
[x, t, U] = varmeledning_theta(D, 1, N, T, M, 0.5, begynnelse_varme, t_vanster, t_hoger, f);
w = U(end,:).';
plot(x, w)
figure;
N = 100; T = 2; dt = 1e-3; M = round(T/dt);
[x, t, U] = varmeledning_theta(D, 1, N, T, M, 0.5, begynnelse_varme, t_vanster, t_hoger, f);
w = U(end,:).';
plot(x, w)
figure;
N = 100; T = 4; dt = 1e-3; M = round(T/dt);
[x, t, U] = varmeledning_theta(D, 1, N, T, M, 0.5, begynnelse_varme, t_vanster, t_hoger, f);
w = U(end,:).';
plot(x, w)
figure;
N = 100; T = 8; dt = 1e-3; M = round(T/dt);
[x, t, U] = varmeledning_theta(D, 1, N, T, M, 0.5, begynnelse_varme, t_vanster, t_hoger, f);
w = U(end,:).';
plot(x, w)

function [x, t, U] = varmeledning_theta(D, L, N, T, M, theta, u0, alpha, beta, f)
    h = L / N;
    x = (0:N)' * h;
    xi = x(2:N);
    dt = T / M;
    t = (0:M)' * dt;

    e = ones(N-1, 1);
    A = D/h^2 * spdiags([e -2*e e], -1:1, N-1, N-1);
    I = speye(N-1);
    P = I - theta * dt * A;
    Q = I + (1-theta) * dt * A; %theta = 1/2 gör att det blir Crank Nicholson

    b = @(t) f(xi,t) + D/h^2 * [alpha(t); zeros(N-3,1); beta(t)]; %f räknar med till och bortföring av värme

    V = u0(xi);
    U = zeros(M+1, N+1);
    U(1,:) = [alpha(0); V; beta(0)]';

    for n = 1:M
        tn = t(n);
        tnp1 = t(n+1);
        V = P \ (Q * V + dt * (theta * b(tnp1) + (1-theta) * b(tn)));
        U(n+1,:) = [alpha(tnp1); V; beta(tnp1)]'; %spara
    end
end