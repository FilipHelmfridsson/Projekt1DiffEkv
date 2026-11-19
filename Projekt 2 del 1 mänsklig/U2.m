clear; close all;


N = 400;
d = 5; gamma = 1; delta = 1;
f    = @(x) d*cos(2*x);
uanalytisk = @(x) (delta - d/12)*cos(4*x) + gamma/4*sin(4*x) + d/12*cos(2*x);

%d
[x,w] = bvp_fdm(N, f, gamma, delta);

w_analytisk = uanalytisk(x);

hold on;
plot(x,w);
plot (x,w_analytisk);
hold off;
figure;

felvektor = abs(w-w_analytisk);

plot(x,felvektor);

%e
N = 50;
fel = []

for i = [1:4]
    N = N*2;
    [x, w] = bvp_fdm(N, f, gamma, delta);
    w_analytisk = uanalytisk(x);
    r = w(2:N) - w_analytisk(2:N);             % felet i de inre punkterna
    fel = [fel, sqrt(sum(abs(r.^2)) / (N-1))];
    
end

noggranhetsordning = log2((fel(1))/fel(2));
noggranhetsordning1 = log2((fel(2))/fel(3));
noggranhetsordning2 = log2((fel(3))/fel(4));

%f)
N = 400;
figure;
g = @(x) exp(-(x-pi).^2);
gamma = 0;
delta = 0;
[x,w] = bvp_fdm(N, g, gamma, delta);
plot(x,w,'r')
hold on;


gamma = 1;
delta = 0;
[x,w] = bvp_fdm(N, g, gamma, delta);
plot(x,w,'b');
hold off;

%Efterssom gamma anger startlutningen som begynnelsevärde, påverkar det hela resten av funktionen.
%Det gör att den ena uppträder betydligt svängigare än den andra.
