function [x, w] = bvp_fdm(N, f, gamma, delta)
% Löser u'' + 16u = f(x), 0 <= x <= 2*pi, u'(0) = gamma, u(2*pi) = delta
% med finita differenser. Returnerar alla N+1 punkter x och lösningen w.

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
