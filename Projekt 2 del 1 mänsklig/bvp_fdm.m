function [x, w] = bvp_fdm(N, f, gamma, delta)

h = 2*pi / N;

x = (0:N)' * h;       %alla punkter
xi = x(2:N);          %inre punkter
n  = N - 1;           %antal obekanta = storleken på ekvationssystemet

A = zeros(n, n);

for j = 1:n

    %diagonalelementet (koefficient framför w_j)
    A(j, j) = -2/h^2 + 16;

    %elementet till vänster (koefficient framför w_(j-1)),
    %finns inte på första raden
    if j > 1
        A(j, j-1) = 1/h^2;
    end

    %elementet till höger (koefficient framför w_{j+1}),
    %finns inte på sista raden
    if j < n
        A(j, j+1) = 1/h^2;
    end

end


A(1, 1) = -2/(3*h^2) + 16;
A(1, 2) =  2/(3*h^2);

q = f(xi);

%Första raden: termen med gamma från steg 4 hamnar i högerledet.
q(1) = q(1) + 2*gamma / (3*h);

q(n) = q(n) - delta / h^2;

wi = A \ q; %lösningen i de inre punkterna w_1 till  w_(N-1)

w0 = (4*wi(1) - wi(2) - 2*h*gamma) / 3;

w = [w0; wi; delta];  %alla N+1 värden: w_0 till w_N

end
