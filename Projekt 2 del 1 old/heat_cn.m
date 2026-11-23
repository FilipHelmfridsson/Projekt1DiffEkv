function [x, w] = heat_cn(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn)
% Löser u_t = D u_xx + f(x,t) på [0,1] x [0,T] med Dirichletrandvillkor
% med Crank–Nicolson. Returnerar x (N+1 punkter) och lösningen w vid t = T.

dx = 1/N;
dt = T/M;
x  = (0:N)' * dx;
xi = x(2:N);
e  = ones(N-1, 1);

A = D/dx^2 * spdiags([e -2*e e], -1:1, N-1, N-1);
I = speye(N-1);
L = I - dt/2*A;
R = I + dt/2*A;

% Källterm + randbidrag: b(t) = f(x_i,t) + D/dx^2 * [TL(t); 0; ...; 0; TR(t)]
b = @(t) ffcn(xi, t) + D/dx^2 * [TLfcn(t); zeros(N-3,1); TRfcn(t)];

wi = uinitfcn(xi);
for n = 0:M-1
    t  = n*dt;
    wi = L \ (R*wi + dt/2*(b(t) + b(t + dt)));
end
w = [TLfcn(T); wi; TRfcn(T)];
end
