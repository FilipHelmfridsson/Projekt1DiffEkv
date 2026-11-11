function [x, w] = bvp_fdm(N, f, gamma, delta)
% Löser u'' + 16u = f(x), 0 <= x <= 2*pi, u'(0) = gamma, u(2*pi) = delta
% med finita differenser. Returnerar alla N+1 punkter x och lösningen w.

h  = 2*pi/N;
x  = (0:N)' * h;
xi = x(2:N);                 % inre punkter x_1 ... x_{N-1}
e  = ones(N-1, 1);

% (w_{j+1} - 2w_j + w_{j-1})/h^2 + 16 w_j = f_j
A = spdiags([e/h^2, (-2/h^2 + 16)*e, e/h^2], -1:1, N-1, N-1);
A(1,1) = -2/(3*h^2) + 16;    % w_0 eliminerad med (6)
A(1,2) =  2/(3*h^2);

q = f(xi);
q(1)   = q(1)   + 2*gamma/(3*h);
q(end) = q(end) - delta/h^2;

wi = A \ q;
w0 = (4*wi(1) - wi(2) - 2*h*gamma) / 3;
w  = [w0; wi; delta];
end
