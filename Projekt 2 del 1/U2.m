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
Ns  = [50 100 200 400 800];
fel = zeros(size(Ns));
for k = 1:length(Ns)
    [x, w] = bvp_fdm(Ns(k), f, gamma, delta);
    r = w(2:end-1) - uana(x(2:end-1));      % inre punkter
    fel(k) = sqrt(mean(r.^2));              % diskret 2-norm (11)
end
ordning = [NaN, log2(fel(1:end-1) ./ fel(2:end))];
table(Ns', fel', ordning', 'VariableNames', {'N', 'Fel', 'Ordning'})

%% f) f = exp(-(x-pi)^2)
g = @(x) exp(-(x - pi).^2);
[x, w1] = bvp_fdm(400, g, 0, 0);
[~, w2] = bvp_fdm(400, g, 1, 0);
figure
plot(x, w1, x, w2)
legend('\gamma = 0, \delta = 0', '\gamma = 1, \delta = 0'), xlabel('x'), ylabel('u')
